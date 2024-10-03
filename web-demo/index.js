
import { threads } from "wasm-feature-detect";
import init, {
  initThreadPool,
  parallel_sum,
  witness_preprocess,
  BfvEncrypt1024,
  BfvEncrypt2048,
  BfvEncrypt4096,
  BfvEncrypt8192,
  BfvEncrypt16384,
  parse_args,
} from "./pkg/web_demo.js";

function setButtonsDisabledState(buttonIds, state) {
  for (let id of buttonIds) {
    let btn = document.getElementById(id);
    if (btn) {
      btn.disabled = state;
    }
  }
}

async function setup() {
  let supportsThreads = await threads();
  if (!supportsThreads) {
    console.error("This browser does not support threads");
    return;
  }
  await init();
  await initThreadPool(2)



  const demoNames = [
    "parallelSum",
    "witnessPreprocess",
    "bfvProve"
  ];

  function setupBtn(id, fn) {
    // Handlers are named in the same way as buttons.
    // let fn = demos[id];

    let button = document.getElementById(id);
    if (button === null) {
      console.error("button with id: ", id, "not found");
      return null;
    }

    // Assign onclick handler + enable the button.
    Object.assign(button, {
      onclick: async () => {
        document.getElementById("loader").hidden = false;
        document.getElementById("testSuccess").checked = false;
        setButtonsDisabledState(demoNames, true);

        console.log("Running: ", id);
        try {
          let results = await fn();
          document.getElementById("testSuccess").checked = true;
          if (results !== undefined) {
            document.getElementById("benchmarkResults").value =
              JSON.stringify(results);
          }
        } catch (error) {
          console.error(`Test Failed: ${error}`);
          document.getElementById("testSuccess").checked = false;
        }
        document.getElementById("loader").hidden = true;
        setButtonsDisabledState(demoNames, false);
      },
      disabled: false,
    });

    return button;
  }


  setupBtn("parallelSum", parallelSum);
  setupBtn("witnessPreprocess", witnessPreprocess);
  setupBtn("bfvProve", bfvProve)

}

setup().then(() => {
  console.log("Setup done")
});


async function parallelSum() {
  const dataSize = 1e7;
  const data = new Float64Array(dataSize);

  // Fill the array with random data
  for (let i = 0; i < dataSize; i++) {
    data[i] = Math.random();
  }
  console.time("parallel_sum");
  // Perform the parallel sum
  const result = parallel_sum(data);
  console.timeEnd("parallel_sum");

}

async function witnessPreprocess() {
  const test_data = [
      14
    // 10, 11, 12, 13, 14
  ];

  for (const i of test_data) {
    const pow2 = 1 << i;
    console.time(`witness_preprocess ${pow2}`);
    const result = witness_preprocess(i);
    console.timeEnd(`witness_preprocess ${pow2}`);
  }

}

async function bfvProve() {
  const test_data = [
    [BfvEncrypt1024.new(), "sk_enc_1024_1x27_65537"],
    [BfvEncrypt2048.new(), "sk_enc_2048_1x52_65537"],
    [BfvEncrypt4096.new(), "sk_enc_4096_2x55_65537"],
    [BfvEncrypt8192.new(), "sk_enc_8192_4x55_65537"],
    [BfvEncrypt16384.new(), "sk_enc_16384_8x54_65537"],
  ];

  for (const [prover, file] of test_data) {
    const response = await fetch(`${file}.json`);
    const jsonData = await response.json();
    // console.log(jsonData);
    let args = parse_args(jsonData);
    // console.log(args);

    console.time(`bfv_gkr setup ${file}`);
    const pk = prover.setup();
    console.timeEnd(`bfv_gkr setup ${file}`);


    console.time(`bfv_gkr prove ${file}`);
    const result = prover.prove(args, pk);
    console.timeEnd(`bfv_gkr prove ${file}`);

    console.log("----------------");
  }


}