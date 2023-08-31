// DOM Elements
let canvas = document.getElementById('canvas');
let ctx = canvas.getContext('2d');
let timeDisplay = document.getElementById('timeDisplay');
let playPauseButton = document.getElementById('playPauseButton');
let restartButton = document.getElementById('restartButton');
let colorScaleCanvas = document.getElementById('colorScale');
let colorScaleCtx = colorScaleCanvas.getContext('2d');

// Simulation variables
let playing = false;
let timeElapsed = 0;
let grid;

// Initialize
function init() {
  // Initialize grid
  grid = new Array(canvas.width);
  for (let i = 0; i < canvas.width; i++) {
    grid[i] = new Array(canvas.height).fill(parseFloat(document.getElementById('T0Rest').value));
  }
  
  // Create an initial hot circle in the center
  let circleRadius = canvas.width * parseFloat(document.getElementById('circleRadius').value) / 10;
  let centerX = canvas.width / 2;
  let centerY = canvas.height / 2;
  
  for (let x = 0; x < canvas.width; x++) {
    for (let y = 0; y < canvas.height; y++) {
      let dx = x - centerX;
      let dy = y - centerY;
      if (dx * dx + dy * dy <= circleRadius * circleRadius) {
        grid[x][y] = parseFloat(document.getElementById('T0Circle').value);
      }
    }
  }
  
  // Draw initial state and color scale
  drawState();
  drawColorScaleBar();
  
  // Attach event listeners
  playPauseButton.addEventListener('click', togglePlayPause);
  restartButton.addEventListener('click', restartSimulation);
}

// Update the simulation
function update() {
  if (playing) {
    // Finite-difference for the heat equation
    let newGrid = JSON.parse(JSON.stringify(grid)); // Deep copy
    let alpha = parseFloat(document.getElementById('thermalDiffusivity').value);
    let dt = parseFloat(document.getElementById('t_stepSize').value);
    let dx = 1;
    
    for (let i = 1; i < grid.length - 1; i++) {
      for (let j = 1; j < grid[0].length - 1; j++) {
        let d2Tdx2 = grid[i+1][j] - 2 * grid[i][j] + grid[i-1][j];
        let d2Tdy2 = grid[i][j+1] - 2 * grid[i][j] + grid[i][j-1];
        newGrid[i][j] = grid[i][j] + alpha * dt * (d2Tdx2 + d2Tdy2) / (dx * dx);
      }
    }
    
    grid = newGrid;
    
    // Update time and display
    timeElapsed += dt;
    timeDisplay.innerHTML = `t = ${timeElapsed.toFixed(2)} s`;
    
    // Draw updated state
    drawState();
  }
  
  // Loop
  requestAnimationFrame(update);
}

// Draw the current state
function drawState() {
  for (let x = 0; x < grid.length; x++) {
    for (let y = 0; y < grid[0].length; y++) {
      let temp = grid[x][y];
      let color = `rgba(${temp}, 0, ${255 - temp}, 1)`;
      ctx.fillStyle = color;
      ctx.fillRect(x, y, 1, 1);
    }
  }
}

// Draw color scale
function drawColorScaleBar() {
  let grad = colorScaleCtx.createLinearGradient(0, 0, 0, colorScaleCanvas.height);
  grad.addColorStop(0, 'blue');
  grad.addColorStop(1, 'red');
  colorScaleCtx.fillStyle = grad;
  colorScaleCtx.fillRect(0, 0, colorScaleCanvas.width, colorScaleCanvas.height);
}

// Toggle play/pause
function togglePlayPause() {
  playing = !playing;
  playPauseButton.innerHTML = playing ? '<i class="fas fa-pause"></i>' : '<i class="fas fa-play"></i>';
}

// Restart simulation
function restartSimulation() {
  init();
  drawState();
  timeElapsed = 0;
  timeDisplay.innerHTML = 't = 0.00 s';
  playing = false;
  playPauseButton.innerHTML = '<i class="fas fa-play"></i>';
}

// Initialize and start
init();
update();