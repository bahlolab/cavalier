// scripts.js

// 1) Read out all panels
function getVariantData() {
  const variantPanels = Array.from(document.querySelectorAll('.variant-info-panel'));
  const genePanels    = Array.from(document.querySelectorAll('.gene-info-panel'));
  const jigvPanels    = Array.from(document.querySelectorAll('.jigv-panel'));

  return variantPanels.map(vp => {
    const id         = vp.dataset.variantId;
    const table_html = vp.innerHTML;
    const gp = genePanels.find(gp => gp.dataset.variantId === id);
    const plot_html  = gp ? gp.innerHTML : '';
    const jp = jigvPanels.find(jp => jp.dataset.variantId === id);
    const jp_html = jp ? jp.innerHTML : '';
    return { id, title: id, table_html, plot_html, jp_html};
  });
}

// 2) Single selectVariant: swaps panels AND jumps IGV
function selectVariant(v) {
  // 1) swap your HTML panels
  document.getElementById("variant-body").innerHTML = v.table_html;
  document.getElementById("gene-body")   .innerHTML = v.plot_html;
  document.getElementById("jigv-body")   .innerHTML = v.jp_html;
}

// 3) Build the sidebar links
function buildSidebar() {
  const list = document.getElementById("variant-list");
  list.innerHTML = "";

  const seen = new Set();
  getVariantData().forEach(v => {
    if (seen.has(v.id)) return;
    seen.add(v.id);
    const a = document.createElement("a");
    a.href        = "#";
    a.textContent = v.title;
    a.addEventListener("click", e => {
      e.preventDefault();
      selectVariant(v);
      // highlight
      document.querySelectorAll("#variant-list a")
              .forEach(x => x.classList.remove("active"));
      a.classList.add("active");
    });
    const li = document.createElement("li");
    li.appendChild(a);
    list.appendChild(li);
  });

  // auto-select the first link
  const first = list.querySelector("a");
  if (first) {
    first.classList.add("active");
    first.click();
  }
}

// 4) Wire it up
document.addEventListener("DOMContentLoaded", buildSidebar);
