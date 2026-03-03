/**
 * InterPro domain visualization - reusable D3 module.
 * Renders protein domains from InterProScan matches.
 * Inspired by GeneNoteBook domain visualization.
 */
(function () {
  'use strict';

  const DEFAULT_MARGIN = { top: 20, right: 20, bottom: 30, left: 100 };
  const ROW_HEIGHT = 30;
  const DEFAULT_WIDTH = 800;

  /**
   * Render InterPro domains on a protein length scale.
   * @param {string} containerSelector - CSS selector for the SVG container (e.g. '#alignment')
   * @param {Array} matches - Array of { name, locations: [{ start, end, evalue }] }
   * @param {number} proteinLength - Total protein length for the x-axis
   * @param {Object} options - Optional { margin, width }
   */
  function renderInterProDomains(containerSelector, matches, proteinLength, options = {}) {
    if (!matches || matches.length === 0) return;

    const margin = options.margin || DEFAULT_MARGIN;
    const width = (options.width || DEFAULT_WIDTH) - margin.left - margin.right;
    const height = matches.length * ROW_HEIGHT;

    const container = d3.select(containerSelector);
    if (container.empty()) return;

    container.selectAll('*').remove();

    container
      .attr('width', width + margin.left + margin.right)
      .attr('height', height + margin.top + margin.bottom);

    const svg = container
      .append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    const xScale = d3.scaleLinear()
      .domain([0, proteinLength])
      .range([0, width]);

    const yScale = d3.scaleBand()
      .domain(matches.map((d) => d.name))
      .range([0, height])
      .padding(0.1);

    svg
      .append('g')
      .attr('transform', `translate(0,${height})`)
      .call(d3.axisBottom(xScale));

    svg.append('g').call(d3.axisLeft(yScale));

    matches.forEach((match) => {
      (match.locations || []).forEach((loc) => {
        const rect = svg
          .append('rect')
          .attr('x', xScale(loc.start))
          .attr('y', yScale(match.name))
          .attr('width', Math.max(1, xScale(loc.end) - xScale(loc.start)))
          .attr('height', yScale.bandwidth())
          .attr('fill', '#69b3a2')
          .attr('opacity', 0.8);

        const tooltip = [
          match.name,
          `${loc.start}–${loc.end}`,
          loc.evalue != null ? `E-value: ${loc.evalue}` : null,
        ].filter(Boolean).join('\n');
        rect.append('title').text(tooltip);
      });
    });
  }

  if (typeof window !== 'undefined') {
    window.renderInterProDomains = renderInterProDomains;
  }
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { renderInterProDomains };
  }
})();
