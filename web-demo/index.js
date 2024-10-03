
import { threads } from "wasm-feature-detect";
import init, {
  initThreadPool,
  prove_encrypt_test
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
  // await initThreadPool(10)



  const demoNames = [
    // "parallelSum",
    "proveEncrypt",
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


  setupBtn("bfvProve", bfvProve)

}

setup().then(() => {
  console.log("Setup done")
});


async function bfvProve() {
  console.time(`total`);
  prove_encrypt_test();
  console.timeEnd(`total`);
}
