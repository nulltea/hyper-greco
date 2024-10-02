build_web:
	wasm-pack build --release --target=web
	find pkg/snippets -type f -iname workerHelpers.worker.js -exec sed -i "s|from '..\/..\/..\/';|from '..\/..\/..\/web_demo.js';|" {} \;

run_web:
	cp ../bfv-gkr/src/data/goldilocks/* dist/
	npm install
	npm run build && serve --config ../serve.json dist/
