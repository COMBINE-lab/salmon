// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// GitHub Pages: site + base. The repo deploys to
// https://combine-lab.github.io/salmon/.
//
// GitHub-Flavored Markdown (tables, strikethrough, autolinks) is enabled by
// default by Astro's Sätteri Markdown processor, and `@astrojs/mdx` inherits
// that setting, so tables render in both `.md` and `.mdx` pages without any
// extra config. (The legacy `markdown.remarkPlugins: [remarkGfm]` API was
// removed here: it now requires the separate `@astrojs/markdown-remark`
// package and is deprecated in favor of configuring the processor directly.)
export default defineConfig({
  site: 'https://combine-lab.github.io',
  base: '/salmon',
  integrations: [
    starlight({
      title: 'salmon',
      description:
        'salmon 2.0 — fast, accurate transcript quantification from RNA-seq reads (Rust).',
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/COMBINE-lab/salmon',
        },
      ],
      editLink: {
        baseUrl: 'https://github.com/COMBINE-lab/salmon/edit/main/website/',
      },
      sidebar: [
        {
          label: 'Getting started',
          items: [
            { label: 'Introduction', slug: 'getting-started/introduction' },
            { label: 'Installation', slug: 'getting-started/installation' },
            { label: 'Quick start', slug: 'getting-started/quick-start' },
          ],
        },
        {
          label: 'Guides',
          items: [
            { label: 'Library types', slug: 'guides/library-types' },
            {
              label: 'Selective alignment & sketch mode',
              slug: 'guides/mapping-modes',
            },
            {
              label: 'Inferential replicates',
              slug: 'guides/inferential-replicates',
            },
            {
              label: 'RAD I/O & deterministic quantification',
              slug: 'guides/rad-and-determinism',
            },
            {
              label: 'Genome-alignment quantification',
              slug: 'guides/genome-projection',
            },
          ],
        },
        {
          label: 'Migrating from C++ salmon',
          items: [{ label: 'What changed in 2.0', slug: 'migrating/from-cpp' }],
        },
        {
          label: 'Reference',
          items: [
            { label: 'Command-line interface', slug: 'reference/cli' },
            { label: 'Output format specification', slug: 'reference/output-formats' },
          ],
        },
      ],
    }),
  ],
});
