function runSimulation() {
  // Show loading indicator
  document.getElementById('loadingIndicator').style.visibility = 'visible';

  setTimeout(() => {
      // Constants (for simplification, I'm hardcoding some values. Ideally, these should be user inputs.)
      const height = parseFloat(document.getElementById("height").value);
      // ... retrieve other inputs ...

      const Nx = 50;
      const Ny = 50;
      const dx = height / (Nx - 1);
      const dy = height / (Ny - 1);
      const abs_coeff = 0.01 + (loading * 2300);

      let q = computePowerSourceDistribution(abs_coeff, beam_radius_m, dx, height, Q, Nx, Ny, dy);
      
      // Visualize power density
      Plotly.newPlot('powerDensityPlot', [{
          z: q,
          type: 'heatmap',
          colorscale: 'Hot'
      }]);
      
      // For this demo, we won't perform the full temperature calculations.
      // Instead, we'll visualize the dummy power distribution.
      Plotly.newPlot('temperaturePlot', [{
          z: q,
          type: 'heatmap',
          colorscale: 'Hot'
      }]);

      // Hide loading indicator after simulations
      document.getElementById('loadingIndicator').style.visibility = 'hidden';
  }, 1000);  // Simulate a 1-second delay for the simulation
}

function computePowerSourceDistribution(abs_coeff, beam_radius_m, dx, height, Q, Nx, Ny, dy) {
  let q = [];
  for (let i = 0; i < Ny; i++) {
      let row = [];
      for (let j = 0; j < Nx; j++) {
          if (Math.sqrt(i*i + j*j) <= beam_radius_m) {
              row.push(Q * Math.exp(-abs_coeff * i * dy));
          } else {
              row.push(0);
          }
      }
      q.push(row);
  }
  return q;
}

// Removed the window.onload function to prevent auto-running the simulation
