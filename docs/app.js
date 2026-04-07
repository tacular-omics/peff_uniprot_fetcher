// Main-thread glue: wires the form to a Web Worker that runs Pyodide.
// The worker does all the heavy lifting (Python + UniProt fetches) so the
// UI stays responsive.

const logEl = document.getElementById("log");
const form = document.getElementById("form");
const goBtn = document.getElementById("go");
const organismInput = document.getElementById("organism");
const reviewedInput = document.getElementById("reviewed");
const downloadSlot = document.getElementById("download-slot");

function log(msg) {
  logEl.textContent += "\n" + msg;
  logEl.scrollTop = logEl.scrollHeight;
}

function clearLog() {
  logEl.textContent = "";
}

const worker = new Worker("./worker.js");

worker.onmessage = (e) => {
  const { type, payload } = e.data;
  if (type === "log") {
    log(payload);
  } else if (type === "ready") {
    log("Pyodide ready.");
    goBtn.disabled = false;
  } else if (type === "result") {
    const { organismId, peff } = payload;
    const blob = new Blob([peff], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    downloadSlot.innerHTML = "";
    const a = document.createElement("a");
    a.href = url;
    a.download = `organism_${organismId}.peff`;
    a.className = "download";
    a.textContent = `Download organism_${organismId}.peff (${(blob.size / 1024).toFixed(1)} KB)`;
    downloadSlot.appendChild(a);
    log(`Done. ${blob.size} bytes.`);
    goBtn.disabled = false;
  } else if (type === "error") {
    log("ERROR: " + payload);
    goBtn.disabled = false;
  }
};

goBtn.disabled = true;
clearLog();
log("Booting Pyodide (first load downloads ~10 MB)...");
worker.postMessage({ type: "init" });

form.addEventListener("submit", (e) => {
  e.preventDefault();
  const organismId = organismInput.value.trim();
  const reviewed = reviewedInput.checked;
  if (!organismId) return;
  goBtn.disabled = true;
  downloadSlot.innerHTML = "";
  log(`\n>>> Generating PEFF for taxon ${organismId} (reviewed=${reviewed})...`);
  worker.postMessage({
    type: "generate",
    payload: { organismId, reviewed },
  });
});
