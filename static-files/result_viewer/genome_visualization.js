// Global variables
let svg, starship_vis, s, scale_vertical, zoom, starship_xAxis, yAxis;
let width = window.innerWidth * 0.95;
let height = config.height;
let chart = document.querySelector('.chart');
let previoust = [], previousb = []; // For strand positioning

// Configuration
const config = {
  margin: { top: 20, right: 20, bottom: 30, left: 50 },
  height: 500,
  defaultOpacity: 0.9,
  transitionDuration: 200
};

// State management
let state = {
  currentScale: null,
  selectedFeature: null,
  isZooming: false
};

// Initialize visualization
function initializeVisualization() {
  setupSvg();
  setupScales();
  setupAxes();
  setupZoom();
  makeResponsive();
  addAccessibility();
  drawFeatures();
}

// Setup functions
function setupSvg() {
  svg = d3.select(".chart")
    .attr("width", width)
    .attr("height", height);
    
  starship_vis = svg.append("g")
    .attr("id", "starship_vis")
    .attr("transform", `translate(${config.margin.left}, ${config.margin.top})`);
}

function setupScales() {
  s = d3.scaleLinear()
    .domain([0, STARSHIP_LENGTH])
    .range([0, width - config.margin.left - config.margin.right]);
    
  scale_vertical = d3.scaleLinear()
    .domain([2, 0])
    .range([0, 50]);
}

// Feature rendering
function drawFeatures() {
  const features = starship_vis.selectAll(".feature")
    .data(feature_data.features)
    .join(
      enter => enterFeatures(enter),
      update => updateFeatures(update),
      exit => exitFeatures(exit)
    );
    
  enhanceFeatureInteractivity(features);
}

function enterFeatures(enter) {
  return enter.append("a")
    .attr("xlink:href", d => d.href)
    .append("rect")
    .attr("class", "feature")
    .call(setFeatureAttributes);
}

function updateFeatures(update) {
  return update.call(setFeatureAttributes);
}

function exitFeatures(exit) {
  return exit.remove();
}

// Helper functions
function setFeatureAttributes(selection) {
  return selection
    .attr("x", d => getFeatureX(d))
    .attr("y", d => getFeatureY(d))
    .attr("width", d => getFeatureWidth(d))
    .attr("height", 5)
    .style("fill", d => getFeatureColor(d))
    .style("opacity", config.defaultOpacity);
}

// Event handlers
function handleFeatureClick(event, d) {
  state.selectedFeature = d;
  showFeatureDetail(d);
}

// Utility functions
function getFeatureColor(d) {
  if (feature_data.feature_id === d.id) return "fuchsia";
  if (d.flag === "N/A" || d.flag === "UNANNOTATED") return "grey";
  // ... rest of color logic ...
}

// Initialize
$(document).ready(() => {
  initializeVisualization();
});

function loadFeatureData() {
  try {
    const featureData = JSON.parse(
      document.getElementById("feature_data").textContent
    );
    return featureData;
  } catch (error) {
    console.error("Error loading feature data:", error);
    return { features: [] };
  }
}

// Add debouncing for performance
const debouncedZoom = _.debounce((event) => {
  const newScale = event.transform.rescaleX(s);
  updateVisualization(newScale);
}, 100);

// Add feature filtering for large datasets
function filterVisibleFeatures(scale) {
  const [start, end] = scale.domain();
  return feature_data.features.filter(d => 
    (d.start >= start && d.start <= end) ||
    (d.stop >= start && d.stop <= end)
  );
}

function makeResponsive() {
  const resizeObserver = new ResizeObserver(_.debounce(() => {
    const newWidth = chart.getBoundingClientRect().width * 0.95;
    const newHeight = chart.getBoundingClientRect().height;
    
    // Update SVG dimensions
    svg.attr("width", newWidth)
       .attr("height", newHeight);
       
    // Update scales
    s.range([0, newWidth - config.margin.left - config.margin.right]);
    
    // Redraw
    redrawVisualization();
  }, 250));
  
  resizeObserver.observe(chart);
  
  // Cleanup
  return () => resizeObserver.disconnect();
}

function addAccessibility() {
  svg.attr("role", "img")
     .attr("aria-label", "Starship visualization")
     .attr("tabindex", 0);
     
  // Add keyboard navigation
  svg.on("keydown", (event) => {
    switch(event.key) {
      case "ArrowRight":
        panRight();
        break;
      case "ArrowLeft":
        panLeft();
        break;
      case "+":
        zoomIn();
        break;
      case "-":
        zoomOut();
        break;
    }
  });
}

function panRight() {
  const currentTransform = d3.zoomTransform(svg.node());
  svg.transition()
    .duration(config.transitionDuration)
    .call(zoom.translateBy, -50, 0);
}

function panLeft() {
  const currentTransform = d3.zoomTransform(svg.node());
  svg.transition()
    .duration(config.transitionDuration)
    .call(zoom.translateBy, 50, 0);
}

function zoomIn() {
  svg.transition()
    .duration(config.transitionDuration)
    .call(zoom.scaleBy, 1.2);
}

function zoomOut() {
  svg.transition()
    .duration(config.transitionDuration)
    .call(zoom.scaleBy, 0.8);
}

function getFeatureX(d) {
  return d.strand === "+" ? s(d.start) : s(d.stop);
}

function getFeatureY(d) {
  return d.strand === "+" ? 
    bottom_strand(previousb, d.start, d.stop) : 
    top_strand(previoust, d.start, d.stop);
}

function getFeatureWidth(d) {
  const start = s(d.start);
  const stop = s(d.stop);
  return Math.abs(stop - start);
}

function showFeatureDetail(feature) {
  // Clear any existing detail view
  d3.select("#feature-detail").remove();
  
  // Create detail panel
  const detailPanel = d3.select("body")
    .append("div")
    .attr("id", "feature-detail")
    .attr("class", "feature-detail-panel");
    
  // Add feature information
  detailPanel.html(`
    <h3>Feature Details</h3>
    <p>ID: ${feature.id}</p>
    <p>Position: ${feature.start}-${feature.stop}</p>
    <p>Strand: ${feature.strand}</p>
    <p>Annotation: ${feature.annotation}</p>
    <p>Flag: ${feature.flag}</p>
    ${feature.public_note ? `<p>Public Note: ${feature.public_note}</p>` : ''}
    ${feature.private_note ? `<p>Private Note: ${feature.private_note}</p>` : ''}
  `);
}

function redrawVisualization() {
  // Update axes
  starship_vis.select(".x-axis").call(starship_xAxis);
  starship_vis.select(".y-axis").call(yAxis);
  
  // Update features
  updateFeatures(s);
}

function updateVisualization(newScale) {
  // Update state
  state.currentScale = newScale;
  
  // Filter features based on visible range
  const visibleFeatures = filterVisibleFeatures(newScale);
  
  // Update visualization
  starship_vis.select(".x-axis").call(starship_xAxis);
  
  // Update features with new positions
  const features = starship_vis.selectAll(".feature")
    .data(visibleFeatures, d => d.id)
    .join(
      enter => enterFeatures(enter),
      update => updateFeatures(update),
      exit => exitFeatures(exit)
    );
}

function setupAxes() {
  // Create X axis
  starship_xAxis = d3.axisBottom(s)
    .ticks(10)
    .tickFormat(d3.format(",.0f"));
    
  // Create Y axis
  yAxis = d3.axisLeft(scale_vertical)
    .ticks(3)
    .tickFormat(d => d === 1 ? "+" : d === 0 ? "-" : "");
    
  // Add axes to visualization
  starship_vis.append("g")
    .attr("class", "x-axis")
    .attr("transform", `translate(0, ${height - config.margin.bottom})`)
    .call(starship_xAxis);
    
  starship_vis.append("g")
    .attr("class", "y-axis")
    .attr("transform", `translate(${-config.margin.left}, 0)`)
    .call(yAxis);
}

function setupZoom() {
  zoom = d3.zoom()
    .scaleExtent([0.5, 20])
    .on("zoom", (event) => {
      if (state.isZooming) return;
      state.isZooming = true;
      
      debouncedZoom(event);
      
      state.isZooming = false;
    });
    
  svg.call(zoom);
  
  // Add zoom controls
  const zoomControls = svg.append("g")
    .attr("class", "zoom-controls")
    .attr("transform", `translate(${width - 60}, 20)`);
    
  zoomControls.append("rect")
    .attr("class", "zoom-reset")
    .attr("width", 40)
    .attr("height", 20)
    .attr("rx", 5)
    .on("click", () => {
      svg.transition()
        .duration(750)
        .call(zoom.transform, d3.zoomIdentity);
    });
}
