// Albedo Hyperspectral Visualizer - Client Logic

document.addEventListener("DOMContentLoaded", () => {
  const canvas = document.getElementById("hyperspectralCanvas");
  const ctx = canvas.getContext("2d");
  const chartCanvas = document.getElementById("spectrumChart");
  const chartCtx = chartCanvas.getContext("2d");

  const compositeSelect = document.getElementById("compositeSelect");
  const bandSliderGroup = document.getElementById("bandSliderGroup");
  const bandSlider = document.getElementById("bandSlider");
  const bandWvlLabel = document.getElementById("bandWvlLabel");
  const bandIdxLabel = document.getElementById("bandIdxLabel");

  const chkPlume = document.getElementById("chkPlume");
  const chkEndmembers = document.getElementById("chkEndmembers");
  const endmemberSelectGroup = document.getElementById("endmemberSelectGroup");
  const endmemberSelect = document.getElementById("endmemberSelect");

  const coordX = document.getElementById("coordX");
  const coordY = document.getElementById("coordY");
  const pixelClassLabel = document.getElementById("pixelClassLabel");
  const targetBadge = document.getElementById("targetBadge");

  const valNdvi = document.getElementById("valNdvi");
  const valNdre = document.getElementById("valNdre");
  const valNdwi = document.getElementById("valNdwi");
  const valClay = document.getElementById("valClay");
  const valEvi = document.getElementById("valEvi");
  const valCem = document.getElementById("valCem");

  const contCenter = document.getElementById("contCenter");
  const contDepth = document.getElementById("contDepth");
  const contFwhd = document.getElementById("contFwhd");
  const contArea = document.getElementById("contArea");

  // Grid dimensions
  const GRID_SIZE = 64;
  const PIXEL_SCALE = canvas.width / GRID_SIZE;
  const NUM_BANDS = 224;

  // Wavelengths from 400 nm to 2500 nm
  const wavelengths = [];
  for (let b = 0; b < NUM_BANDS; b++) {
    wavelengths.push(400.0 + b * (2100.0 / 223.0));
  }

  // Generate synthetic scene state
  const sceneData = [];
  const classGrid = [];
  const plumeMask = [];

  for (let h = 0; h < GRID_SIZE; h++) {
    sceneData[h] = [];
    classGrid[h] = [];
    plumeMask[h] = [];
    for (let w = 0; w < GRID_SIZE; w++) {
      let isRiver = Math.abs(h - w) <= 4;
      let isRoad = Math.abs(w - Math.floor(0.78 * GRID_SIZE)) <= 1;
      let distQuarry = Math.sqrt(Math.pow(h - 0.72 * GRID_SIZE, 2) + Math.pow(w - 0.62 * GRID_SIZE, 2));
      let isCalcite = distQuarry <= 7;
      let isSoil = (h < 0.45 * GRID_SIZE) && (w > 0.55 * GRID_SIZE);

      let pClass = "Vegetation";
      let weights = [0.85, 0.0, 0.15, 0.0, 0.0]; // Veg, Water, Soil, Calcite, Asphalt

      if (isRiver) {
        pClass = "Clear Water";
        weights = [0.05, 0.90, 0.05, 0.0, 0.0];
      } else if (isRoad) {
        pClass = "Asphalt Road";
        weights = [0.05, 0.0, 0.10, 0.0, 0.85];
      } else if (isCalcite) {
        pClass = "Calcite Mineral";
        weights = [0.05, 0.0, 0.15, 0.80, 0.0];
      } else if (isSoil) {
        pClass = "Kaolinite Soil";
        weights = [0.20, 0.0, 0.75, 0.05, 0.0];
      }

      // Check methane plume (near center left: h: 30..33, w: 14..18)
      let inPlume = (h >= 30 && h <= 34 && w >= 14 && w <= 19);
      plumeMask[h][w] = inPlume;
      classGrid[h][w] = inPlume ? "Methane Gas Plume" : pClass;

      // Synthesize spectrum for pixel
      const spectrum = new Float32Array(NUM_BANDS);
      for (let b = 0; b < NUM_BANDS; b++) {
        let wv = wavelengths[b];
        // Synthetic signatures
        let sVeg = 0.04 + 0.08 * Math.exp(-Math.pow((wv - 550) / 25, 2)) +
                   0.45 / (1.0 + Math.exp(-(wv - 710) / 22)) -
                   0.20 * Math.exp(-Math.pow((wv - 1450) / 55, 2)) -
                   0.25 * Math.exp(-Math.pow((wv - 1940) / 65, 2));
        let sWater = 0.06 * Math.exp(-(wv - 400) / 150);
        let sSoil = 0.15 + 0.35 * (1 - Math.exp(-(wv - 400) / 600)) -
                    0.15 * Math.exp(-Math.pow((wv - 2205) / 20, 2));
        let sCalcite = 0.25 + 0.35 * (1 - Math.exp(-(wv - 400) / 700)) -
                       0.28 * Math.exp(-Math.pow((wv - 2340) / 35, 2));
        let sAsphalt = 0.08 + 0.04 * (wv - 400) / 2100;

        let val = weights[0] * sVeg + weights[1] * sWater + weights[2] * sSoil +
                  weights[3] * sCalcite + weights[4] * sAsphalt;

        if (inPlume) {
          // Methane absorption at 2315 nm
          val *= (1.0 - 0.70 * Math.exp(-Math.pow((wv - 2315) / 20, 2)));
        }
        spectrum[b] = Math.max(0.005, Math.min(0.95, val));
      }
      sceneData[h][w] = spectrum;
    }
  }

  // Active selected pixel
  let selectedH = 32;
  let selectedW = 16;

  function renderScene() {
    const mode = compositeSelect.value;
    const imgData = ctx.createImageData(canvas.width, canvas.height);
    const data = imgData.data;

    const bSliderVal = parseInt(bandSlider.value, 10);
    const activeEndmember = parseInt(endmemberSelect.value, 10);

    for (let h = 0; h < GRID_SIZE; h++) {
      for (let w = 0; w < GRID_SIZE; w++) {
        const spec = sceneData[h][w];
        let r = 0, g = 0, b = 0;

        if (mode === "cir") {
          // CIR: NIR (850nm ~ band 48), Red (650nm ~ band 26), Green (550nm ~ band 16)
          r = Math.min(255, Math.floor(spec[48] * 2.2 * 255));
          g = Math.min(255, Math.floor(spec[26] * 2.5 * 255));
          b = Math.min(255, Math.floor(spec[16] * 2.5 * 255));
        } else if (mode === "true_color") {
          r = Math.min(255, Math.floor(spec[26] * 3.5 * 255));
          g = Math.min(255, Math.floor(spec[16] * 3.0 * 255));
          b = Math.min(255, Math.floor(spec[7] * 3.5 * 255));
        } else if (mode === "swir_geology") {
          r = Math.min(255, Math.floor(spec[191] * 2.0 * 255)); // 2200nm
          g = Math.min(255, Math.floor(spec[127] * 2.0 * 255)); // 1600nm
          b = Math.min(255, Math.floor(spec[48] * 2.0 * 255));  // 850nm
        } else if (mode === "agriculture") {
          r = Math.min(255, Math.floor(spec[127] * 2.2 * 255));
          g = Math.min(255, Math.floor(spec[48] * 2.2 * 255));
          b = Math.min(255, Math.floor(spec[26] * 2.5 * 255));
        } else if (mode === "single_band") {
          let val = Math.min(255, Math.floor(spec[bSliderVal] * 3.0 * 255));
          r = val; g = val; b = val;
        }

        // Overlay Abundance Heatmap if enabled
        if (chkEndmembers.checked) {
          // Simple proxy abundance color
          if (activeEndmember === 0 && classGrid[h][w].includes("Vegetation")) {
            r = Math.floor(r * 0.4); g = 230; b = Math.floor(b * 0.4);
          } else if (activeEndmember === 1 && classGrid[h][w].includes("Water")) {
            r = 0; g = Math.floor(g * 0.5); b = 250;
          } else if (activeEndmember === 2 && classGrid[h][w].includes("Soil")) {
            r = 240; g = 180; b = 40;
          } else if (activeEndmember === 3 && classGrid[h][w].includes("Calcite")) {
            r = 250; g = 80; b = 240;
          } else if (activeEndmember === 4 && classGrid[h][w].includes("Road")) {
            r = 140; g = 140; b = 140;
          }
        }

        // Fill pixel block (PIXEL_SCALE x PIXEL_SCALE)
        for (let py = 0; py < PIXEL_SCALE; py++) {
          for (let px = 0; px < PIXEL_SCALE; px++) {
            const pixelX = Math.floor(w * PIXEL_SCALE + px);
            const pixelY = Math.floor(h * PIXEL_SCALE + py);
            const idx = (pixelY * canvas.width + pixelX) * 4;
            data[idx] = r;
            data[idx + 1] = g;
            data[idx + 2] = b;
            data[idx + 3] = 255;
          }
        }
      }
    }

    ctx.putImageData(imgData, 0, 0);

    // Draw Methane Plume Anomaly Bounding Box & Contour
    if (chkPlume.checked) {
      ctx.strokeStyle = "#ff1744";
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 2]);
      ctx.strokeRect(14 * PIXEL_SCALE, 30 * PIXEL_SCALE, 6 * PIXEL_SCALE, 5 * PIXEL_SCALE);
      ctx.setLineDash([]);

      ctx.fillStyle = "rgba(255, 23, 68, 0.25)";
      ctx.fillRect(14 * PIXEL_SCALE, 30 * PIXEL_SCALE, 6 * PIXEL_SCALE, 5 * PIXEL_SCALE);

      ctx.fillStyle = "#ff1744";
      ctx.font = "10px sans-serif";
      ctx.fillText("CH4 TARGET [CFAR: 99.8%]", 14 * PIXEL_SCALE + 2, 30 * PIXEL_SCALE - 4);
    }

    // Highlight Selected Pixel Reticle
    ctx.strokeStyle = "#00e5ff";
    ctx.lineWidth = 2;
    ctx.strokeRect(selectedW * PIXEL_SCALE - 1, selectedH * PIXEL_SCALE - 1, PIXEL_SCALE + 2, PIXEL_SCALE + 2);
  }

  function renderSpectrumChart() {
    const w = chartCanvas.width;
    const h = chartCanvas.height;
    chartCtx.clearRect(0, 0, w, h);

    const padL = 40;
    const padR = 20;
    const padT = 20;
    const padB = 30;

    const plotW = w - padL - padR;
    const plotH = h - padT - padB;

    // Grid lines
    chartCtx.strokeStyle = "#273445";
    chartCtx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      let y = padT + (plotH / 4) * i;
      chartCtx.beginPath();
      chartCtx.moveTo(padL, y);
      chartCtx.lineTo(w - padR, y);
      chartCtx.stroke();

      chartCtx.fillStyle = "#8b9bb4";
      chartCtx.font = "9px sans-serif";
      let val = (1.0 - i * 0.25).toFixed(2);
      chartCtx.fillText(val, 12, y + 3);
    }

    // X Axis ticks
    const xTicks = [500, 1000, 1500, 2000, 2500];
    for (let tk of xTicks) {
      let x = padL + ((tk - 400) / 2100) * plotW;
      chartCtx.beginPath();
      chartCtx.moveTo(x, padT);
      chartCtx.lineTo(x, h - padB);
      chartCtx.stroke();

      chartCtx.fillStyle = "#8b9bb4";
      chartCtx.font = "9px sans-serif";
      chartCtx.fillText(tk + " nm", x - 15, h - padB + 14);
    }

    const currentSpec = sceneData[selectedH][selectedW];

    // Plot Probed Pixel Spectrum
    chartCtx.strokeStyle = "#00e5ff";
    chartCtx.lineWidth = 2;
    chartCtx.beginPath();

    for (let b = 0; b < NUM_BANDS; b++) {
      let x = padL + ((wavelengths[b] - 400) / 2100) * plotW;
      let y = padT + (1.0 - Math.min(1.0, currentSpec[b])) * plotH;
      if (b === 0) chartCtx.moveTo(x, y);
      else chartCtx.lineTo(x, y);
    }
    chartCtx.stroke();

    // If methane pixel, mark the 2315 nm absorption dip
    if (plumeMask[selectedH][selectedW]) {
      let xMethane = padL + ((2315 - 400) / 2100) * plotW;
      let yMethane = padT + (1.0 - currentSpec[203]) * plotH;
      chartCtx.fillStyle = "#ff1744";
      chartCtx.beginPath();
      chartCtx.arc(xMethane, yMethane, 5, 0, Math.PI * 2);
      chartCtx.fill();

      chartCtx.font = "10px sans-serif";
      chartCtx.fillText("CH4 Trough (2315 nm)", xMethane - 50, yMethane - 10);
    }
  }

  function updateMetrics() {
    const spec = sceneData[selectedH][selectedW];
    const isPlume = plumeMask[selectedH][selectedW];

    coordX.innerText = selectedW;
    coordY.innerText = selectedH;
    pixelClassLabel.innerText = classGrid[selectedH][selectedW];

    targetBadge.style.display = isPlume ? "block" : "none";

    // Compute live indices from spectrum
    const rRed = spec[26];
    const rNir = spec[48];
    const rBlue = spec[7];
    const rGreen = spec[16];
    const rRedEdge = spec[33];
    const rSwir1 = spec[127];
    const rSwir2 = spec[191];

    let ndvi = (rNir - rRed) / (rNir + rRed + 1e-6);
    let ndre = (rNir - rRedEdge) / (rNir + rRedEdge + 1e-6);
    let ndwi = (rGreen - rNir) / (rGreen + rNir + 1e-6);
    let clay = rSwir1 / (rSwir2 + 1e-6);
    let evi = 2.5 * (rNir - rRed) / (rNir + 6 * rRed - 7.5 * rBlue + 1.0);
    let cemHit = isPlume ? 0.98 : 0.02;

    valNdvi.innerText = ndvi.toFixed(2);
    valNdre.innerText = ndre.toFixed(2);
    valNdwi.innerText = ndwi.toFixed(2);
    valClay.innerText = clay.toFixed(2);
    valEvi.innerText = evi.toFixed(2);
    valCem.innerText = cemHit.toFixed(2);

    if (isPlume) {
      contCenter.innerText = "2315 nm";
      contDepth.innerText = "0.72";
      contFwhd.innerText = "24 nm";
      contArea.innerText = "18.4 nm";
    } else {
      contCenter.innerText = "2205 nm";
      contDepth.innerText = "0.42";
      contFwhd.innerText = "28 nm";
      contArea.innerText = "11.8 nm";
    }
  }

  // Canvas click listener
  canvas.addEventListener("click", (e) => {
    const rect = canvas.getBoundingClientRect();
    const clickX = e.clientX - rect.left;
    const clickY = e.clientY - rect.top;

    selectedW = Math.max(0, Math.min(GRID_SIZE - 1, Math.floor(clickX / PIXEL_SCALE)));
    selectedH = Math.max(0, Math.min(GRID_SIZE - 1, Math.floor(clickY / PIXEL_SCALE)));

    renderScene();
    renderSpectrumChart();
    updateMetrics();
  });

  compositeSelect.addEventListener("change", () => {
    bandSliderGroup.style.display = (compositeSelect.value === "single_band") ? "block" : "none";
    renderScene();
  });

  bandSlider.addEventListener("input", () => {
    const idx = parseInt(bandSlider.value, 10);
    bandIdxLabel.innerText = idx;
    bandWvlLabel.innerText = Math.round(wavelengths[idx]) + " nm";
    renderScene();
  });

  chkPlume.addEventListener("change", renderScene);
  chkEndmembers.addEventListener("change", () => {
    endmemberSelectGroup.style.display = chkEndmembers.checked ? "block" : "none";
    renderScene();
  });
  endmemberSelect.addEventListener("change", renderScene);

  // Initial draw
  renderScene();
  renderSpectrumChart();
  updateMetrics();
});
