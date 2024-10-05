build_web:
	RUSTFLAGS='-C target-feature=+atomics,+bulk-memory,+mutable-globals' rustup run nightly-2024-09-06 wasm-pack build --target web -- -Z build-std=panic_abort,std;
	find pkg/snippets -type f -iname workerHelpers.worker.js -exec sed -i "s|from '..\/..\/..\/';|from '..\/..\/..\/web_demo.js';|" {} \;

run_web:
	cp ../bfv-gkr/src/data/goldilocks/* dist/
	npm install
	npm run build && serve --config ../serve.json dist/
