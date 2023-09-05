document.addEventListener("DOMContentLoaded", function() {
  const canvas = document.getElementById('canvas');
  const ctx = canvas.getContext('2d');
  const size = 300;
  const dx = 1;
  let running = false;
  let timeElapsed = 0;

  let T, circleRadius, centerX, centerY, dt, alpha;

  function initializeSimulation() {
    const T0Circle = parseFloat(document.getElementById("T0Circle").value);
    const T0Rest = parseFloat(document.getElementById("T0Rest").value);
    circleRadius = Math.floor(size * parseFloat(document.getElementById("circleRadius").value));
    dt = parseFloat(document.getElementById("t_stepSize").value);
    alpha = parseFloat(document.getElementById("thermalDiffusivity").value);

    if (alpha < 0) {
      alert("Thermal diffusivity cannot be negative.");
      return;
    }

    T = Array.from({ length: size }, () => Array(size).fill(T0Rest));
    centerX = Math.floor(size / 2);
    centerY = Math.floor(size / 2);

    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        const distance = Math.sqrt(Math.pow(i - centerX, 2) + Math.pow(j - centerY, 2));
        if (distance <= circleRadius) {
          T[i][j] = T0Circle;
        }
      }
    }
  }

  const colors = [
    '#a50026',
    '#d73027',
    '#f46d43',
    '#fdae61',
    '#fee090',
    '#e0f3f8',
    '#abd9e9',
    '#74add1',
    '#4575b4',
    '#313695'
  ];

  function populateScaleBar() {
    const scaleWrapper = document.getElementById('scaleWrapper');
    scaleWrapper.innerHTML = "<span id='scaleMax'>1000°C</span>";

    for (let i = colors.length - 1; i >= 0; i--) {
      const scaleBlock = document.createElement('div');
      scaleBlock.classList.add('scaleBlock');
      scaleBlock.style.height = `${Math.floor(300 / colors.length)}px`;
      scaleBlock.style.backgroundColor = colors[i];
      scaleWrapper.appendChild(scaleBlock);
    }

    scaleWrapper.innerHTML += "<span id='scaleMin'>0°C</span>";
  }

  function getColor(temp) {
    const index = Math.floor(temp / 1000 * (colors.length - 1));
    return colors[Math.min(colors.length - 1, Math.max(0, index))];
  }

  function step() {
    let newT = JSON.parse(JSON.stringify(T));

    for (let i = 1; i < size - 1; i++) {
      for (let j = 1; j < size - 1; j++) {
        let d2Tdx2 = (T[i+1][j] - 2 * T[i][j] + T[i-1][j]) / (dx * dx);
        let d2Tdy2 = (T[i][j+1] - 2 * T[i][j] + T[i][j-1]) / (dx * dx);
        newT[i][j] = T[i][j] + alpha * dt * (d2Tdx2 + d2Tdy2);
      }
    }

    T = newT;

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        const color = getColor(T[i][j]);
        ctx.fillStyle = color;
        ctx.fillRect(i, j, 1, 1);
      }
    }

    timeElapsed += dt;
    document.getElementById("timeElapsed").textContent = timeElapsed.toFixed(1);
  }

  document.getElementById('restartButton').addEventListener('click', function() {
    running = false;
    document.getElementById('playPauseButton').innerHTML = '<i class="fas fa-play"></i>';
    initializeSimulation();
    timeElapsed = 0;  // Reset the elapsed time
    document.getElementById("timeElapsed").textContent = timeElapsed.toFixed(1);
    step();
  });

  document.getElementById('playPauseButton').addEventListener('click', function() {
    running = !running;
    this.innerHTML = running ? '<i class="fas fa-pause"></i>' : '<i class="fas fa-play"></i>';
  });

  function mainLoop() {
    if (running) {
      step();
    }
    requestAnimationFrame(mainLoop);
  }

  initializeSimulation();
  populateScaleBar();
  mainLoop();
});
