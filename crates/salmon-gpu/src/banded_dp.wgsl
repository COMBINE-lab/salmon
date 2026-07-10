// Banded affine-gap DP, one alignment per invocation, band-relative storage.
//
// A mirror of `reference::run_dp`: same band, same recurrence, same `mqe`-or-max
// extraction, in integer arithmetic, so the result is deterministic and
// bit-identical to the CPU reference.
//
// Each invocation keeps only the *band* of each DP row (width 2w+1), not the
// whole target row, in four small private arrays sized at MAX bandwidth. This
// keeps per-invocation memory tiny (good GPU occupancy) and bounds the inner
// loop to the band (not the whole target). The host routes any task whose
// bandwidth exceeds MAX_W to the CPU reference; there is no target-length cap.

const NEG_INF: i32 = -0x40000000;
const HALF_NEG: i32 = -0x20000000; // NEG_INF / 2, the "is -inf" guard
const WILDCARD: u32 = 4u;
const MAX_W: i32 = 40;
const BW: u32 = 82u; // 2 * MAX_W + 2, the max band width stored per row

struct Params {
    match_score: i32,
    mismatch_pen: i32,
    gap_open_pen: i32,
    gap_extend_pen: i32,
    bandwidth: i32,
    n_tasks: u32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> tasks: array<vec4<u32>>;
@group(0) @binding(2) var<storage, read> qbuf: array<u32>; // one DNA5 base per element
@group(0) @binding(3) var<storage, read> tbuf: array<u32>;
@group(0) @binding(4) var<storage, read_write> scores: array<i32>;

fn sat_sub(a: i32, b: i32) -> i32 {
    if (a <= HALF_NEG) {
        return NEG_INF;
    }
    return a - b;
}

fn sub_score(qb: u32, tb: u32) -> i32 {
    if (qb == WILDCARD || tb == WILDCARD) {
        return -params.gap_extend_pen;
    }
    if (qb == tb) {
        return params.match_score;
    }
    return -params.mismatch_pen;
}

// Band rows, indexed relative to each row's left band column.
var<private> h_prev: array<i32, BW>; // H(iq-1, .) relative to lo_prev
var<private> h_cur: array<i32, BW>; // H(iq, .)   relative to lo
var<private> f_prev: array<i32, BW>; // F(iq-1, .) relative to lo_prev
var<private> f_cur: array<i32, BW>; // F(iq, .)   relative to lo

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let task = gid.x;
    if (task >= params.n_tasks) {
        return;
    }
    let tmeta = tasks[task];
    let q_off = tmeta.x;
    let q_len = tmeta.y;
    let t_off = tmeta.z;
    let t_len = tmeta.w;

    let go = params.gap_open_pen;
    let ge = params.gap_extend_pen;
    let w = params.bandwidth;
    if (q_len == 0u || t_len == 0u || w > MAX_W) {
        scores[task] = NEG_INF; // host recomputes these on the CPU reference
        return;
    }
    let qlen = i32(q_len);
    let tlen = i32(t_len);

    var global_max: i32 = NEG_INF;
    var mqe: i32 = NEG_INF;
    var lo_prev: i32 = 0;
    var hi_prev: i32 = -1; // row -1 is virtual; handled analytically below

    for (var iq: u32 = 0u; iq < q_len; iq = iq + 1u) {
        let iqi = i32(iq);
        let lo = max(0, iqi - w);
        let hi = min(tlen - 1, iqi + w);
        let qb = qbuf[q_off + iq];

        // H(iq, it-1) running across the band, seeded at the band's left edge.
        var h_left: i32 = select(NEG_INF, -(go + (iqi + 1) * ge), lo == 0);
        var e_run: i32 = NEG_INF;

        for (var it: i32 = lo; it <= hi; it = it + 1) {
            let rel = u32(it - lo);

            // Exact per-cell band membership (matches the reference).
            let r = iqi + it;
            var st: i32 = 0;
            var en: i32 = tlen - 1;
            if (st < r - qlen + 1) { st = r - qlen + 1; }
            if (en > r) { en = r; }
            let blo = (r - w + 1) >> 1u;
            let bhi = (r + w) >> 1u;
            if (st < blo) { st = blo; }
            if (en > bhi) { en = bhi; }
            if (!(it >= st && it <= en)) {
                h_cur[rel] = NEG_INF;
                f_cur[rel] = NEG_INF;
                e_run = NEG_INF;
                h_left = NEG_INF;
                continue;
            }

            // Diagonal source H(iq-1, it-1).
            var diag: i32;
            if (it == 0) {
                diag = select(-(go + iqi * ge), 0, iq == 0u); // H(iq-1, -1)
            } else if (iq == 0u) {
                diag = -(go + it * ge); // H(-1, it-1)
            } else {
                let j = it - 1;
                if (j >= lo_prev && j <= hi_prev) {
                    diag = h_prev[u32(j - lo_prev)];
                } else {
                    diag = NEG_INF;
                }
            }

            let s = sub_score(qb, tbuf[t_off + u32(it)]);
            var m: i32 = NEG_INF;
            if (diag > HALF_NEG) {
                m = diag + s;
            }

            // E: gap in query (consume target), from the left neighbour.
            e_run = max(sat_sub(h_left, go + ge), sat_sub(e_run, ge));

            // F: gap in target (consume query), from the cell above.
            var h_up: i32;
            var f_up: i32;
            if (iq == 0u) {
                h_up = -(go + (it + 1) * ge); // H(-1, it)
                f_up = NEG_INF;
            } else if (it >= lo_prev && it <= hi_prev) {
                h_up = h_prev[u32(it - lo_prev)];
                f_up = f_prev[u32(it - lo_prev)];
            } else {
                h_up = NEG_INF;
                f_up = NEG_INF;
            }
            let f = max(sat_sub(h_up, go + ge), sat_sub(f_up, ge));
            f_cur[rel] = f;

            var h = max(max(m, e_run), f);
            if (h < NEG_INF) { h = NEG_INF; }
            h_cur[rel] = h;
            h_left = h;

            if (h > global_max) { global_max = h; }
            if (iq + 1u == q_len && h > mqe) { mqe = h; }
        }

        // Roll the current band into the previous-row band for the next row.
        let width = u32(hi - lo + 1);
        for (var k: u32 = 0u; k < width; k = k + 1u) {
            h_prev[k] = h_cur[k];
            f_prev[k] = f_cur[k];
        }
        lo_prev = lo;
        hi_prev = hi;
    }

    var result: i32;
    if (mqe > NEG_INF) {
        result = mqe;
    } else {
        result = max(global_max, 0);
    }
    scores[task] = result;
}
