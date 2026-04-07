// Web Worker that hosts Pyodide + peff_uniprot_fetcher.
//
// Running in a worker has two key benefits:
//   1. The main thread stays responsive during multi-minute runs.
//   2. Synchronous XMLHttpRequest is permitted in workers, which is what
//      peff_uniprot_fetcher/_client.py uses under Pyodide so that the
//      existing synchronous pipeline keeps working unchanged.

importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js");

let pyodide = null;
let ready = false;

function logToMain(msg) {
  self.postMessage({ type: "log", payload: String(msg) });
}

async function init() {
  if (ready) {
    self.postMessage({ type: "ready" });
    return;
  }
  logToMain("Loading Pyodide runtime...");
  pyodide = await loadPyodide({
    indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/",
    stdout: logToMain,
    stderr: logToMain,
  });

  logToMain("Installing micropip...");
  await pyodide.loadPackage("micropip");

  logToMain("Installing peff_uniprot_fetcher and dependencies (from PyPI)...");
  await pyodide.runPythonAsync(`
import micropip
# Install runtime deps from PyPI. All are pure-Python wheels.
await micropip.install([
    "httpx",
    "pefftacular",
    "psimodpy",
    "unimodpy",
    "uniprotptmpy",
])
# Install our local wheel (served alongside this worker).
await micropip.install("./peff_uniprot_fetcher-0.1.0-py3-none-any.whl")

import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
from peff_uniprot_fetcher._web import generate_peff_string
`);

  ready = true;
  self.postMessage({ type: "ready" });
}

async function generate(organismId, reviewed) {
  if (!ready) {
    self.postMessage({ type: "error", payload: "Pyodide not ready yet" });
    return;
  }
  try {
    pyodide.globals.set("_organism_id", String(organismId));
    pyodide.globals.set("_reviewed", Boolean(reviewed));
    const peff = await pyodide.runPythonAsync(`
generate_peff_string(_organism_id, reviewed=_reviewed)
`);
    self.postMessage({
      type: "result",
      payload: { organismId, peff },
    });
  } catch (err) {
    self.postMessage({ type: "error", payload: err && err.message ? err.message : String(err) });
  }
}

self.onmessage = async (e) => {
  const { type, payload } = e.data;
  if (type === "init") {
    try {
      await init();
    } catch (err) {
      self.postMessage({ type: "error", payload: "Init failed: " + (err && err.message ? err.message : String(err)) });
    }
  } else if (type === "generate") {
    await generate(payload.organismId, payload.reviewed);
  }
};
