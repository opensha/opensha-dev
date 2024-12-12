package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.CubedGriddedRegion;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.LineString;
import org.opensha.commons.geo.json.Geometry.GeometryCollection;
import org.opensha.commons.geo.json.Geometry.Polygon;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_FaultCubeAssociations;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.erf.ETAS.SeisDepthDistribution;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class GridCubeRatePlot {

	public static void main(String[] args) throws IOException {
		double lonSpacingKM = LocationUtils.horzDistance(new Location(0d, 0d), new Location(0d, 1d));
		double gridSpacing = 0.1;
		double cellSpacingKM = gridSpacing * lonSpacingKM;
		double upDepth = 0d;
		double lowDepth = 15d;
//		int cellsBefore = 2;
//		int cellsAfter = 2;
//		double shadeInCellFract = 0.3;
		int cellsBefore = 2;
		int cellsAfter = 2;
		double shadeInCellFract = 0.0;
//		int cellsBefore = 1;
//		int cellsAfter = 2;
//		double shadeInCellFract = 0.0;
		boolean vertical = false;
		double dip, rake;
		Location traceStart, traceEnd;
		if (vertical) {
			dip = 90d;
			rake = 0d;
			traceStart = new Location(0d, 0d);
			traceEnd = new Location(1d, 0d);
		} else {
			// set the dip such that the association 2 cells away is zero
			double assocWidth = NSHM23_SingleRegionGridSourceProvider.DEFAULT_MAX_FAULT_NUCL_DIST;
			double fromSideOfCell = assocWidth % cellSpacingKM;
			System.out.println("Calculated dist from site of cell: "+(float)fromSideOfCell);
			Preconditions.checkState(fromSideOfCell*2 < cellSpacingKM);
			double horzWidth = cellSpacingKM - 2*fromSideOfCell;
			double height = lowDepth - upDepth;
			dip = Math.toDegrees(Math.atan(height/horzWidth));
			double traceLon = -0.5*gridSpacing + fromSideOfCell/lonSpacingKM;
			traceStart = new Location(0d, traceLon);
			traceEnd = new Location(1d, traceLon);
			rake = 0d;
		}
		
		File outputDir = new File(FIGURES_DIR, "crustal_grid");
		
		Region sectReg = new Region(new Location(traceStart.lat, traceStart.lon-1d), new Location(traceEnd.lat, traceEnd.lon+1d));
		GriddedRegion gridReg = new GriddedRegion(sectReg, gridSpacing, GriddedRegion.ANCHOR_0_0);
		
		GutenbergRichterMagFreqDist supraMFD = new GutenbergRichterMagFreqDist(1d, 0.1, 6.55, 7.55, 11);
		EvenlyDiscretizedFunc gridRefMFD = FaultSysTools.initEmptyMFD(5.05, 7.55);
		GutenbergRichterMagFreqDist gridMFD = new GutenbergRichterMagFreqDist(1d, 1d, gridRefMFD.getMinX(), gridRefMFD.getMaxX(), gridRefMFD.size());
		System.out.println("Have "+gridReg.getNodeCount()+" total grid nodes");
		double rateEach = 1e-2;
		double sumRate = rateEach*gridReg.getNodeCount();
		gridMFD.scaleToCumRate(5.05, sumRate);
		
		FeatureProperties props = new FeatureProperties();
		props.set(GeoJSONFaultSection.FAULT_ID, 0);
		props.set(GeoJSONFaultSection.FAULT_NAME, "Test Fault");
		props.set(GeoJSONFaultSection.DIP, dip);
		props.set(GeoJSONFaultSection.RAKE, rake);
		props.set(GeoJSONFaultSection.LOW_DEPTH, lowDepth);
		props.set(GeoJSONFaultSection.UPPER_DEPTH, upDepth);
		props.set(GeoJSONFaultSection.SLIP_RATE, 10d);
		LineString geom = new LineString(traceStart, traceEnd);
//		props.set(GeoJSONFaultSection.PROXY, true);
//		Geometry geom = new LineString(traceStart, traceEnd);
//		Polygon proxyZone = new Polygon(new Region(traceStart, new Location(traceEnd.lat, traceEnd.lon+gridSpacing)));
//		geom = new GeometryCollection(geom, proxyZone);
		Feature sectFeature = new Feature(geom, props);
		GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(sectFeature);
		
		Location centerLoc = new Location(0.5*(sectReg.getMinLat()+sectReg.getMaxLat()),
				0.5*(sectReg.getMinLon()+sectReg.getMaxLon()));
		int centerGridIndex = gridReg.indexForLocation(centerLoc);
		
		List<Integer> plotCells = new ArrayList<>();
		for (int i=0; i<cellsBefore; i++)
			plotCells.add(0, gridReg.indexForLocation(new Location(centerLoc.lat, centerLoc.lon-gridSpacing*(i+1d))));
		plotCells.add(centerGridIndex);
		for (int i=0; i<cellsAfter; i++)
			plotCells.add(gridReg.indexForLocation(new Location(centerLoc.lat, centerLoc.lon+gridSpacing*(i+1d))));
		System.out.println("Using cylidrical lon spacing of "+(float)lonSpacingKM+" km");
		double plotMinLon = Double.POSITIVE_INFINITY;
		double plotMaxLon = Double.NEGATIVE_INFINITY;
		System.out.println("Will plot across "+plotCells.size()+" cells:");
		List<Range> cellXRanges = new ArrayList<>();
		for (int cell : plotCells) {
			Location loc = gridReg.getLocation(cell);
			System.out.println("\t"+loc);
			double cellStart = loc.lon-0.5*gridSpacing;
			double cellEnd = loc.lon+0.5*gridSpacing;
			plotMinLon = Math.min(plotMinLon, cellStart);
			plotMaxLon = Math.max(plotMaxLon, cellEnd);
			cellXRanges.add(new Range(cellStart*lonSpacingKM, cellEnd*lonSpacingKM));
		}
		
		Range xRange = new Range(plotMinLon*lonSpacingKM + shadeInCellFract*cellSpacingKM, plotMaxLon*lonSpacingKM - shadeInCellFract*cellSpacingKM);
		
		double[] rupMags = new double[supraMFD.size()];
		double[] rupAreas = new double[supraMFD.size()];
		double[] rupLengths = new double[supraMFD.size()];
		double[] rupRates = new double[supraMFD.size()];
		List<List<Integer>> sectsForRup = new ArrayList<>(rupRates.length);
		for (int i=0; i<rupRates.length; i++) {
			rupMags[i] = supraMFD.getX(i);
			rupAreas[i] = sect.getArea(false);
			rupLengths[i] = sect.getTraceLength();
			rupRates[i] = supraMFD.getY(i);
			sectsForRup.add(List.of(0));
		}
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.builder(List.of(sect), sectsForRup)
				.rupMags(rupMags)
				.rupLengths(rupLengths)
				.rupAreas(rupAreas)
				.build();
		FaultSystemSolution sol = new FaultSystemSolution(rupSet, rupRates);
		
		CubedGriddedRegion cgr = new CubedGriddedRegion(gridReg);
		
		NSHM23_FaultCubeAssociations cubeAssoc = new NSHM23_FaultCubeAssociations(rupSet, cgr,
				NSHM23_SingleRegionGridSourceProvider.DEFAULT_MAX_FAULT_NUCL_DIST);
		rupSet.addModule(cubeAssoc);
		
		double[] pdf = new double[gridReg.getNodeCount()];
		double[] fracSS = new double[gridReg.getNodeCount()];
		double[] fracN = new double[gridReg.getNodeCount()];
		double[] fracR = new double[gridReg.getNodeCount()];
		for (int i=0; i<pdf.length; i++) {
			pdf[i] = 1d/pdf.length;
			fracSS[i] = 1d;
		}
		
		boolean preserveMFD = false;
		SeisDepthDistribution seisDepthDistribution = new SeisDepthDistribution();
		double delta=2;
		HistogramFunction binnedDepthDistFunc = new HistogramFunction(1d, 12,delta);
		for(int i=0;i<binnedDepthDistFunc.size();i++) {
			double prob = seisDepthDistribution.getProbBetweenDepths(binnedDepthDistFunc.getX(i)-delta/2d,binnedDepthDistFunc.getX(i)+delta/2d);
			binnedDepthDistFunc.set(i,prob);
		}
		NSHM23_SingleRegionGridSourceProvider gridProv = new NSHM23_SingleRegionGridSourceProvider(
				sol, cubeAssoc, pdf, gridMFD, preserveMFD, binnedDepthDistFunc, fracSS, fracN, fracR, Map.of(TectonicRegionType.ACTIVE_SHALLOW, sectReg));
		
		System.out.println("Plot X Range: "+xRange);
		
		double cubeSpacingX = cgr.getCubeLatLonSpacing()*lonSpacingKM;
		double cubeSpacingZ = cgr.getCubeDepthDiscr();
		List<int[]> cellCubes = new ArrayList<>();
		double minCubeX = Double.POSITIVE_INFINITY;
		double maxCubeX = Double.NEGATIVE_INFINITY;
		double minCubeZ = Double.POSITIVE_INFINITY;
		double maxCubeZ = Double.NEGATIVE_INFINITY;
		for (int cell : plotCells) {
			int[] cubes = cgr.getCubeIndicesForGridCell(cell);
			cellCubes.add(cubes);
			for (int cube : cubes) {
				Location cubeLoc = cgr.getCubeLocationForIndex(cube);
				double x = cubeLoc.lon*lonSpacingKM;
				double z = cubeLoc.depth;
				minCubeX = Math.min(minCubeX, x);
				maxCubeX = Math.max(maxCubeX, x);
				minCubeZ = Math.min(minCubeZ, z);
				maxCubeZ = Math.max(maxCubeZ, z);
			}
		}
		System.out.println("Cube xSpacing: "+(float)cubeSpacingX);
		System.out.println("Cube zSpacing: "+(float)cubeSpacingZ);
		System.out.println("Cube xRange: ["+(float)minCubeX+", "+(float)maxCubeX+"]");
		System.out.println("Cube zRange: ["+(float)minCubeZ+", "+(float)maxCubeZ+"]");
		int nx = (int)((maxCubeX - minCubeX)/cubeSpacingX) + 1;
		int nz = (int)((maxCubeZ - minCubeZ)/cubeSpacingZ) + 1;
		EvenlyDiscrXYZ_DataSet assocXYZ = new EvenlyDiscrXYZ_DataSet(nx, nz, minCubeX, minCubeZ, cubeSpacingX, cubeSpacingZ);
		Preconditions.checkState((float)assocXYZ.getMaxX() == (float)maxCubeX,
				"x mismatch; nx=%s, cubeMaxX=%s, xyzMaxX=%s", nx, maxCubeX, assocXYZ.getMaxX());
		Preconditions.checkState((float)assocXYZ.getMaxY() == (float)maxCubeZ,
				"z mismatch; nz=%s, cubeMaxZ=%s, xyzMaxY=%s", nz, maxCubeZ, assocXYZ.getMaxY());
		List<List<Integer>> xyzMappedCubes = new ArrayList<>();
		for (int i=0; i<assocXYZ.size(); i++)
			xyzMappedCubes.add(new ArrayList<>());
		for (int[] cubes : cellCubes) {
			for (int cube : cubes) {
				Location cubeLoc = cgr.getCubeLocationForIndex(cube);
				double x = cubeLoc.lon*lonSpacingKM;
				double z = cubeLoc.depth;
				int xyzIndex = assocXYZ.indexOf(x, z);
				Preconditions.checkArgument(xyzIndex >= 0);
				xyzMappedCubes.get(xyzIndex).add(cube);
			}
		}
		int minNumMapped = Integer.MAX_VALUE;
		int maxNumMapped = 0;
		for (List<Integer> mapped : xyzMappedCubes) {
			minNumMapped = Integer.min(minNumMapped, mapped.size());
			maxNumMapped = Integer.max(maxNumMapped, mapped.size());
		}
		System.out.println("Cube to xyz mapping range: ["+minNumMapped+", "+maxNumMapped+"]");
		
		XY_DataSet sectXY = new DefaultXY_DataSet();
		Location sectTopLoc = sect.getFaultTrace().first();
		Location sectBotLoc = sect.getFaultSurface(1d).getEvenlyDiscritizedLowerEdge().first();
		sectXY.set(sectTopLoc.lon*lonSpacingKM, sectTopLoc.depth);
		sectXY.set(sectBotLoc.lon*lonSpacingKM, sectBotLoc.depth);
		
		CPT assocCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(0d, 1d);
		assocCPT.setPreferredTickInterval(0.1);
		calcXYZValues(assocXYZ, xyzMappedCubes, new Function<Integer, Double>() {
			
			@Override
			public Double apply(Integer cubeIndex) {
				double[] weights = cubeAssoc.getOrigSectDistWeightsAtCube(cubeIndex);
				if (weights == null)
					return 0d;
				return StatUtils.sum(weights);
			}
		}, false); // average
//		List<Double> cellAssocs = aggregateToCells(assocXYZ, cellXRanges, false);
		List<Double> cellAssocs = new ArrayList<>();
		for (int cell : plotCells)
			cellAssocs.add(cubeAssoc.getNodeFraction(cell));
		System.out.println("Cell associations: "+cellAssocs);
		
		DecimalFormat fractDF = new DecimalFormat("0.00");
		
		XYZPlotSpec assocPlot = buildPlot(assocXYZ, cellXRanges, cellAssocs, new FormatFunc(fractDF),
				sectXY, assocCPT, "Fractional Association", " ");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		Range yRange = new Range(CELL_TOP_Z, assocXYZ.getMaxY()+0.5*cubeSpacingZ);
		gp.drawGraphPanel(assocPlot, false, false, xRange, yRange);
		
		int width = plotCells.size() > 4 ? 1000 : 800;
		
		PlotUtils.writePlots(outputDir, "cube_assoc", gp, width, false, true, true, false);
		
		EvenlyDiscrXYZ_DataSet totalNuclXYZ = assocXYZ.copy();
		EvenlyDiscrXYZ_DataSet m7NuclXYZ = assocXYZ.copy();
		calcXYZValues(totalNuclXYZ, xyzMappedCubes, new Function<Integer, Double>() {
			
			@Override
			public Double apply(Integer cubeIndex) {
//				return gridProv.getTotalMFD_ForCube(cubeIndex).calcSumOfY_Vals();
				return gridProv.getGriddedSeisMFD_ForCube(cubeIndex).calcSumOfY_Vals();
			}
		}, true); // sum them
		calcXYZValues(m7NuclXYZ, xyzMappedCubes, new Function<Integer, Double>() {
			
			@Override
			public Double apply(Integer cubeIndex) {
//				return gridProv.getTotalMFD_ForCube(cubeIndex).getCumRate(gridRefMFD.getClosestXIndex(7.01));
				return gridProv.getGriddedSeisMFD_ForCube(cubeIndex).getCumRate(gridRefMFD.getClosestXIndex(7.01));
			}
		}, true); // sum them
		List<Double> totCellRates = new ArrayList<>();
		List<Double> m7CellRates = new ArrayList<>();
		totalNuclXYZ.log10();
		m7NuclXYZ.log10();
		for (int cell : plotCells) {
			IncrementalMagFreqDist cellMFD = gridProv.getMFD(cell);
			totCellRates.add(cellMFD.calcSumOfY_Vals());
			m7CellRates.add(cellMFD.getCumRate(cellMFD.getClosestXIndex(7.01)));
		}
		System.out.println("Cell total nucl rates: "+totCellRates);
		System.out.println("Cell M>7 nucl rates: "+m7CellRates);
		for (int i=0; i<totCellRates.size(); i++) {
			totCellRates.set(i, Math.log10(totCellRates.get(i)));
			m7CellRates.set(i, Math.log10(m7CellRates.get(i)));
		}
		
		CPT nuclRateCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(-7, -2d);
		nuclRateCPT.setPreferredTickInterval(0.5);
		Function<Double, String> rateFormat = new Function<Double, String>() {
			
			private DecimalFormat eDF = new DecimalFormat("0.000E0");
			
			@Override
			public String apply(Double t) {
				// back to linear
				double linear = Math.pow(10, t);
				return eDF.format(linear).replace("E", "e");
			}
		};
		
		XYZPlotSpec cellRatePlot = buildPlot(totalNuclXYZ, cellXRanges, totCellRates, rateFormat,
				sectXY, nuclRateCPT, "Gridded M>5 Nucleation Rate", " ");
		
		gp.drawGraphPanel(cellRatePlot, false, false, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "cube_nucl_rate_tot", gp, width, false, true, true, false);
		
		XYZPlotSpec cellRateM7Plot = buildPlot(m7NuclXYZ, cellXRanges, m7CellRates, rateFormat,
				sectXY, nuclRateCPT, "Gridded M>7 Nucleation Rate", " ");
		
		gp.drawGraphPanel(cellRateM7Plot, false, false, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "cube_nucl_rate_m7", gp, width, false, true, true, false);
	}
	
	private static void calcXYZValues(EvenlyDiscrXYZ_DataSet xyz, List<List<Integer>> xyzMappedCubes,
			Function<Integer, Double> cubeValFunc, boolean sum) {
		for (int i=0; i<xyz.size(); i++) {
			List<Integer> mappedCubes = xyzMappedCubes.get(i);
			Preconditions.checkState(!mappedCubes.isEmpty());
			double sumVals = 0d;
			for (int cube : mappedCubes)
				sumVals += cubeValFunc.apply(cube);
			if (sum)
				xyz.set(i, sumVals);
			else
				// take average
				xyz.set(i, sumVals/mappedCubes.size());
		}
	}

	private static final double CELL_TOP_Z = -1.8d;
	private static final double CELL_LINE_SHIFT = 0.05;
	private static final double CELL_LINE_TOP_Z = CELL_TOP_Z+CELL_LINE_SHIFT;
	private static final double CELL_LINE_BOT_Z = -CELL_LINE_SHIFT;
	
	private static List<Double> aggregateToCells(EvenlyDiscrXYZ_DataSet xyz, List<Range> cellXRanges, boolean sum) {
		double[] cellSums = new double[cellXRanges.size()];
		int[] cellCounts = new int[cellXRanges.size()];
		for (int i=0; i<xyz.size(); i++) {
			double x = xyz.getPoint(i).getX();
			double val = xyz.get(i);
			int cellIndex = -1;
			for (int j=0; j<cellXRanges.size(); j++) {
				if (cellXRanges.get(j).contains(x)) {
					cellIndex = j;
					break;
				}
			}
			Preconditions.checkState(cellIndex >= 0, "Cube with x=%s not contained in cellXRanges=%s", x, cellXRanges);
			cellSums[cellIndex] += val;
			cellCounts[cellIndex]++;
		}
		if (!sum) {
			// take average
			List<Double> ret = new ArrayList<>();
			for (int i=0; i<cellSums.length; i++) {
				Preconditions.checkState(cellCounts[i] > 0);
				ret.add(cellSums[i]/cellCounts[i]);
			}
			return ret;
		}
		return Doubles.asList(cellSums);
	}
	
	private static class FormatFunc implements Function<Double, String> {
		
		private DecimalFormat df;

		public FormatFunc(DecimalFormat df) {
			this.df = df;
		}

		@Override
		public String apply(Double t) {
			return df.format(t);
		}
		
	}
	
	private static XYZPlotSpec buildPlot(EvenlyDiscrXYZ_DataSet xyz, List<Range> cellXRanges, List<Double> cellValues,
			Function<Double, String> cellDF, XY_DataSet sectXY, CPT cpt, String zLabel, String title) {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		List<XYTextAnnotation> cellAnns = cellDF == null ? null : new ArrayList<>();
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, cellValues.size() > 3 ? 16 : 20);
		
		for (int i=0; i<cellValues.size(); i++) {
			Color color = cpt.getColor(cellValues.get(i).floatValue());
			Range cellXRange = cellXRanges.get(i);
			
			XY_DataSet outline = new DefaultXY_DataSet();
			double lower = cellXRange.getLowerBound();
			double upper = cellXRange.getUpperBound();
			if (i == 0)
				lower += CELL_LINE_SHIFT;
			else if (i == cellValues.size()-1)
				upper -= CELL_LINE_SHIFT;
			outline.set(lower, CELL_LINE_BOT_Z);
			outline.set(upper, CELL_LINE_BOT_Z);
			outline.set(upper, CELL_LINE_TOP_Z);
			outline.set(lower, CELL_LINE_TOP_Z);
			outline.set(outline.get(0));
			
			funcs.add(outline);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			
			funcs.add(0, outline);
			chars.add(0, new PlotCurveCharacterstics(PlotLineType.POLYGON_SOLID, 1f, color));
			
			if (cellDF != null) {
				XYTextAnnotation ann = new XYTextAnnotation(cellDF.apply(cellValues.get(i)), cellXRange.getCentralValue(), 0.5*CELL_TOP_Z);
				ann.setTextAnchor(TextAnchor.CENTER);
				ann.setFont(annFont);;
				ann.setBackgroundPaint(new Color(255, 255, 255, 127));
				cellAnns.add(ann);
			}
		}
		
		funcs.add(sectXY);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 7f, Color.WHITE));
		
		funcs.add(sectXY);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		XYZPlotSpec plot = new XYZPlotSpec(xyz, funcs, chars, cpt, title, "X (km)", "Depth (km)", zLabel);
		plot.setYAxisInverted(true);
		plot.setPlotAnnotations(cellAnns);
		return plot;
	}

}
