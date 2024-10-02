build_web:
	wasm-pack build --release --target=web
	find pkg/snippets -type f -iname workerHelpers.worker.js -exec sed -i "s|from '..\/..\/..\/';|from '..\/..\/..\/web_demo.js';|" {} \;
