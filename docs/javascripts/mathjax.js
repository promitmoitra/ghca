// MathJax 3 config for MkDocs Material + pymdownx.arithmatex (generic mode).
// arithmatex wraps math in \(...\) / \[...\] with class "arithmatex"; MathJax
// typesets those. Re-typeset on Material's instant-navigation page swaps.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
