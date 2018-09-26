/**
* Single-cell equivalence class
**/
#[derive(Debug)]
pub struct SCEQClass {
    // transcripts defining this eq. class
    pub transcripts: Vec<u32>,
    // cell barcodes ids where this eq. class appears
    pub cells: Vec<u32>,
    pub cell_counts: Vec<u32>,
    // umis with multiplicities
    // the k-mer should be a k-mer class eventually
    pub umis: Vec<(String, u32)>,
}

impl SCEQClass {
    pub fn new(
        transcripts: Vec<u32>,
        cells: Vec<u32>,
        cell_counts: Vec<u32>,
        umis: Vec<(String, u32)>,
    ) -> SCEQClass {
        SCEQClass {
            transcripts: transcripts,
            cells: cells,
            cell_counts: cell_counts,
            umis: umis,
        }
    }
}

// processed experiment
#[derive(Debug)]
pub struct ProcessedExp {
    pub transcripts: Vec<String>,
    pub num_cells: u64,
    pub barcodes: Vec<String>,
    pub eq_classes: Vec<SCEQClass>,
}

/**
 * Single-cell equivalence class
 **/
#[derive(Debug)]
pub struct CellEQClass<'a> {
    // transcripts defining this eq. class
    pub transcripts: &'a Vec<u32>,
    // umis with multiplicities
    // the k-mer should be a k-mer class eventually
    pub umis: Vec<(String, u32)>,
}

#[derive(Debug)]
pub struct CellExp<'a> {
    pub transcripts: &'a Vec<String>,
    pub eq_classes: Vec<CellEQClass<'a>>,
}
