import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Pseudovirus deep mutational scanning of CHIKV envelope proteins",
  description:
    "How mutations to the CHIKV envelope proteins affect entry in various cells and Mxra8 binding",
  base: "/CHIKV-181-25-E-DMS/",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/CHIKV-181-25-E-DMS" }],
    footer: {
      message: "Copyright © 2025-present Xiaohui Ju and Jesse Bloom",
    },
  },
});
