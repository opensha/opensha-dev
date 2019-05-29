package scratch.kevin.simulators.plots;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.DAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectDAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SlipLengthScalingPlot extends AbstractPlot {
	
	private RSQSimSubSectionMapper mapper;
	private double minMag;
	private SlipAlongSectAlgorithm[] slipAlgs;
	private Map<SlipAlongSectAlgorithm, DefaultXY_DataSet> slipAlgResultsMap;
	
	private RSQSimEvent exampleEvent;

	public static final double LEN_COMP_DDW = 12;
	public static final int MIN_VALS_FOR_AVERAGE = 20;

	public SlipLengthScalingPlot(RSQSimSubSectionMapper mapper, double minMag, SlipAlongSectAlgorithm... slipAlgs) {
		if (slipAlgs == null || slipAlgs.length == 0)
			slipAlgs = SlipAlongSectAlgorithm.values();
		this.mapper = mapper;
		this.minMag = minMag;
		this.slipAlgs = slipAlgs;
		
		mapper.trackSlipOnSections();
		
		slipAlgResultsMap = new HashMap<>();
		for (SlipAlongSectAlgorithm alg : slipAlgs)
			slipAlgResultsMap.put(alg, new DefaultXY_DataSet());
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		if (e.getMagnitude() < minMag)
			return;
		Preconditions.checkState(e instanceof RSQSimEvent);
		List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings((RSQSimEvent)e);
		
		checkExampleCandidate((RSQSimEvent)e, mappings);
		
		for (SlipAlongSectAlgorithm slipAlg : slipAlgs) {
			double sumLength = 0d;
			double areaWeightedSlip = 0d;
			double sumArea = 0d;
			for (List<SubSectionMapping> bundle : mappings) {
				for (SubSectionMapping mapping : bundle) {
					sumLength += mapping.getLengthForSlip(slipAlg);
					double area = mapping.getAreaForAverageSlip(slipAlg);
					double slip = mapping.getAverageSlip(slipAlg);
					areaWeightedSlip += slip*area;
					sumArea += area;
				}
			}
			double aveSlip = areaWeightedSlip == 0 ? 0 : areaWeightedSlip / sumArea;
			slipAlgResultsMap.get(slipAlg).set(sumLength, aveSlip);
		}
	}
	
	private static boolean EXAMPLE_DEBUG = false;
	
	private void checkExampleCandidate(RSQSimEvent event, List<List<SubSectionMapping>> allMappings) {
		if (exampleEvent != null)
			return;
		if (allMappings.size() != 1)
			// we only want one parent section
			return;
		List<SubSectionMapping> mappings = allMappings.get(0);
		if (mappings.size() < 3)
			// we want at least 3 total mapped subsections
			return;
//		System.out.println("Have "+mappings.size()+" for example");
		for (SimulatorElement elem : event.getAllElements())
			if (elem.getFocalMechanism().getDip() < 89d)
				// we only want strike slip
				return;
//		List<List<SubSectionMapping>> filteredMappings = mapper.getFilteredSubSectionMappings(event, 0.2);
//		if (filteredMappings.size() != 1 || filteredMappings.get(0).size() != mappings.size()-1)
//			// we want 3 fully mapped subsections and 1 partially ruptured
//			return;
		boolean D = EXAMPLE_DEBUG;
		if (D) System.out.println("Initial example candidate with "+mappings.size()+" sections at t="+event.getTimeInYears());
		double totalLen = 0d;
		double slippedLen = 0d;
		double midSlippedLen = 0d;
		double surfSlippedLen = 0d;
		DAS_Record prevSurfDAS = null;
		for (int i=0; i<mappings.size(); i++) {
			boolean internal = i > 0 && i < mappings.size()-1;
			SubSectionMapping mapping = mappings.get(i);
			if (mapping.isReversed())
				return;
			double myTotLen = mapping.getLengthForSlip(SlipAlongSectAlgorithm.MID_SEIS_FULL_LEN);
			double mySlippedLen = mapping.getLengthForSlip(SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN);
			double myMidSlippedLen = mapping.getLengthForSlip(SlipAlongSectAlgorithm.MID_SEIS_MID_SLIPPED_LEN);
			double mySurfSlippedLen = mapping.getLengthForSlip(SlipAlongSectAlgorithm.MID_SEIS_SURF_SLIP_LEN);
			if (internal) {
				// make sure the internal subsections are fully ruptured
				if ((float)myTotLen != (float)mySlippedLen || (float)mySlippedLen != (float)myMidSlippedLen) {
					if (D) System.out.println("internal not fully ruptured");
					return;
				}
			}
			if (mySurfSlippedLen > 0) {
				DAS_Record surfDAS = mapping.getDASforSlip(SlipAlongSectAlgorithm.MID_SEIS_SURF_SLIP_LEN);
				// this one has surface slip, make sure it's contiguous
				if (surfSlippedLen > 0 && prevSurfDAS == null) {
					// we had surface slip earlier, but not in the previous section, and it just started again
					if (D) System.out.println("surf slipped earlier, then stopped, then started");
					return;
				}
				if (prevSurfDAS != null) {
					// we had surface slip in the previous rupture
					// make sure we're starting at zero
					if ((float)surfDAS.startDAS > 0f) {
						if (D) System.out.println("surf slipped earlier, but we start at "+(float)surfDAS.startDAS);
						return;
					}
					// make sure that the previous DAS went to the end of the section
					if ((float)prevSurfDAS.endDAS == (float)totalLen) {
						if (D) System.out.println("surf slipped earlier, but didn't end at end");
						return;
					}
				}
				
				prevSurfDAS = surfDAS;
			} else {
				prevSurfDAS = null;
			}
			totalLen += myTotLen;
			slippedLen += mySlippedLen;
			midSlippedLen += myMidSlippedLen;
			surfSlippedLen += mySurfSlippedLen;
		}
		if (D) System.out.println("Checking lenghts: tot="+(float)totalLen+"\tslipped="+(float)slippedLen
				+"\tmidSlipped="+(float)midSlippedLen+"\tsurfSlipped="+(float)surfSlippedLen);
		if (slippedLen < 5 || midSlippedLen < 5 || surfSlippedLen < 5)
			// make sure we have some slip for all
			return;
		if (slippedLen - midSlippedLen < 3d)
			// at least 2 km more of slipped length than mid slipped length
			return;
		if (totalLen - slippedLen < 3d)
			// at least 2km more of total length than slipped length
			return;
		if (slippedLen - surfSlippedLen < 2d)
			// some but not all surface slip
			return;
		if (midSlippedLen == surfSlippedLen)
			return;
		if (totalLen < 30 || totalLen > 60)
			return;
		// we have a match!
		System.out.println("Found an example event!");
		exampleEvent = event;
		if (D) {
			System.out.println("Plotting now...");
			try {
				plotExample();
			} catch (IOException e) {
				 e.printStackTrace();
			}
			System.out.println("DONE");
		}
	}
	
	private static final int max_scatter_points = 100000;

	@Override
	public void finalizePlot() throws IOException {
		String name = getCatalogName();
		File outputDir = getOutputDir();
		int plotWidth = getPlotWidth();
		int plotHeight = getPlotHeight();
		
		double maxY = 0d;
		double maxX = 0d;
		for (XY_DataSet scatter : slipAlgResultsMap.values()) {
			maxY = Math.max(maxY, scatter.getMaxY());
			maxX = Math.max(maxX, scatter.getMaxX());
		}
		Range xRange = new Range(0, maxX*1.1);
		Range yRange = new Range(0, Math.ceil(maxY));
		
		// comparisons
		List<DiscretizedFunc> compFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> compChars = new ArrayList<>();
		ScalingRelationships[] compScales = {
				ScalingRelationships.SHAW_2009_MOD,
				ScalingRelationships.ELLSWORTH_B,
				ScalingRelationships.HANKS_BAKUN_08,
				ScalingRelationships.ELLB_SQRT_LENGTH,
				ScalingRelationships.SHAW_CONST_STRESS_DROP
		};
		Color[] compColors = {
				Color.BLUE.darker(),
				Color.RED.darker(),
				Color.GREEN.darker(),
				Color.MAGENTA.darker(),
				Color.ORANGE.darker()
		};
		for (int i=0; i<compScales.length; i++) {
			compFuncs.add(new ArbitrarilyDiscretizedFunc(compScales[i].getShortName()));
			compChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, compColors[i]));
		}
		
		double origDDW = LEN_COMP_DDW*1e3;
    	for (double len=1d; len<=xRange.getUpperBound(); len++) {
    		double length = len*1e3;
    		double area = length*LEN_COMP_DDW*1e3;
    		for (int i=0; i<compScales.length; i++)
    			compFuncs.get(i).set(len, compScales[i].getAveSlip(area, length, origDDW));
    	}
		
		for (SlipAlongSectAlgorithm slipAlg : slipAlgs) {
			String prefix = getOutputPrefix()+"_"+slipAlg.name();
			XY_DataSet scatter = slipAlgResultsMap.get(slipAlg);
			
			// build 2D hist
			int nx = 51, ny = 51;
			double gridSpacingX = maxX/(double)nx;
			double xyzMinX = 0.5*gridSpacingX;
			double gridSpacingY = maxY/(double)ny;
			double xyzMinY = 0.5*gridSpacingY;
			
			EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(nx, ny, xyzMinX, xyzMinY, gridSpacingX, gridSpacingY);
			
			for (Point2D pt : scatter) {
				int index = xyz.indexOf(pt.getX(), pt.getY());
				if (index < 0 || index >= xyz.size())
					throw new IllegalStateException("Scatter point not in XYZ range. x: "
								+pt.getX()+" ["+xyz.getMinX()+" "+xyz.getMaxX()
							+"], y: "+pt.getY()+" ["+xyz.getMinY()+" "+xyz.getMaxY()+"]");
				xyz.set(index, xyz.get(index)+1);
			}
			
			// build average func
			DiscretizedFunc averageFunc = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<xyz.getNumX(); i++) {
				double totNum = 0d;
				double avgY = 0d;
				for (int j=0; j<xyz.getNumY(); j++) {
					double y = xyz.getY(j);
					double num = xyz.get(i, j);
					avgY += y*num;
					totNum += num;
				}
				if (totNum < MIN_VALS_FOR_AVERAGE) {
					if (averageFunc.size() == 0)
						// continue until we get enough data
						continue;
					// wer're at the end, stop
					break;
				}
				avgY /= totNum;
				averageFunc.set(xyz.getX(i), avgY);
			}
			averageFunc.setName("Average");
			
			XY_DataSet plotScatter = scatter;
			if (scatter.size() > max_scatter_points) {
				System.out.println("Filtering slip-length scatter from "+scatter.size()+" to ~"+max_scatter_points+" points");
//				plotScatter = new DefaultXY_DataSet();
//				Random r = new Random();
//				for (int i=0; i<max_scatter_points; i++)
//					plotScatter.set(scatter.get(r.nextInt(scatter.size())));
				plotScatter = downsampleByMag(scatter, true, max_scatter_points);
				System.out.println("Filter done (random mag-dependent sample): "+plotScatter.size());
			}
			List<XY_DataSet> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(plotScatter);
			plotScatter.setName(name);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
			
			if (averageFunc.size() > 0) {
				funcs.add(averageFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
			}
			
			funcs.addAll(compFuncs);
			chars.addAll(compChars);
			
			System.out.println("Scatter y range: "+scatter.getMinY()+" "+scatter.getMaxY());
			
//			EvenlyDiscretizedFunc wcFunc = null, elbFunc = null, hbFunc = null;
//			PlotCurveCharacterstics wcChar = null, elbChar = null, hbChar = null;
//			if (!meanSlip) {
//				WC1994_MagAreaRelationship wc = new WC1994_MagAreaRelationship();
//				wcFunc = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), 1000);
//				for (int i=0; i<wcFunc.size(); i++)
//					wcFunc.set(i, wc.getMedianMag(wcFunc.getX(i)));
//				wcFunc.setName("W-C 1994");
//				funcs.add(wcFunc);
//				wcChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED.darker());
//				chars.add(wcChar);
//				
//				Ellsworth_B_WG02_MagAreaRel elb = new Ellsworth_B_WG02_MagAreaRel();
//				elbFunc = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), 1000);
//				for (int i=0; i<elbFunc.size(); i++)
//					elbFunc.set(i, elb.getMedianMag(elbFunc.getX(i)));
//				elbFunc.setName("EllsworthB");
//				funcs.add(elbFunc);
//				elbChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE.darker());
//				chars.add(elbChar);
//				
//				HanksBakun2002_MagAreaRel hb = new HanksBakun2002_MagAreaRel();
//				hbFunc = new EvenlyDiscretizedFunc(scatter.getMinX(), scatter.getMaxX(), 1000);
//				for (int i=0; i<hbFunc.size(); i++)
//					hbFunc.set(i, hb.getMedianMag(hbFunc.getX(i)));
//				hbFunc.setName("H-B 2002");
//				funcs.add(hbFunc);
//				hbChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker());
//				chars.add(hbChar);
//			}
			
			String title = "Slip-Length Scaling";
			String xAxisLabel = slipAlg+" (km)";
			String yAxisLabel = "Mean Mid-Seismogenic Slip (m)";
			
			PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
			plot.setLegendVisible(funcs.size() > 1);
			
			HeadlessGraphPanel gp = buildGraphPanel();
			gp.drawGraphPanel(plot, false, false, xRange, yRange);
			gp.getChartPanel().setSize(plotWidth, plotHeight);
			gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			
			// convert to density
			for (int i=0; i<xyz.size(); i++) {
				// convert to density
				double binWidth = gridSpacingX;
				double binHeight = gridSpacingY;
				double area = binWidth * binHeight;
				xyz.set(i, xyz.get(i)*area);
			}
			xyz.scale(1d/xyz.getSumZ());
			
			// set all zero to NaN so that it will plot white
			for (int i=0; i<xyz.size(); i++) {
				if (xyz.get(i) == 0)
					xyz.set(i, Double.NaN);
			}
			xyz.log10();
			
			double minZ = Double.POSITIVE_INFINITY;
			double maxZ = Double.NEGATIVE_INFINITY;
			for (int i=0; i<xyz.size(); i++) {
				double val = xyz.get(i);
				if (!Doubles.isFinite(val))
					continue;
				if (val < minZ)
					minZ = val;
				if (val > maxZ)
					maxZ = val;
			}
			
			System.out.println("MinZ: "+minZ);
			System.out.println("MaxZ: "+maxZ);
			
			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
			if ((float)minZ == (float)maxZ)
				cpt = cpt.rescale(minZ, minZ*2);
			else if (!Doubles.isFinite(minZ))
				cpt = cpt.rescale(0d, 1d);
			else
				cpt = cpt.rescale(minZ, maxZ);
			cpt.setNanColor(Color.WHITE);
			
			String zAxisLabel = "Log10(Density)";
			
			XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, title, xAxisLabel, yAxisLabel, zAxisLabel);
			// add the mean and comparison funcs, but remove the scatter
			funcs.remove(0);
			chars.remove(0);
			xyzSpec.setXYElems(funcs);
			xyzSpec.setXYChars(chars);
//			if (!meanSlip) {
//				// add W-C
//				funcs = Lists.newArrayList();
//				chars = Lists.newArrayList();
//				ArbitrarilyDiscretizedFunc logWCFunc = new ArbitrarilyDiscretizedFunc();
//				for (Point2D pt : wcFunc)
//					logWCFunc.set(Math.log10(pt.getX()), pt.getY());
//				funcs.add(logWCFunc);
//				chars.add(wcChar);
//				// add EllB
//				ArbitrarilyDiscretizedFunc logEllBFunc = new ArbitrarilyDiscretizedFunc();
//				for (Point2D pt : elbFunc)
//					logEllBFunc.set(Math.log10(pt.getX()), pt.getY());
//				funcs.add(logEllBFunc);
//				chars.add(elbChar);
//				// add HB
//				ArbitrarilyDiscretizedFunc logHBFunc = new ArbitrarilyDiscretizedFunc();
//				for (Point2D pt : hbFunc)
//					logHBFunc.set(Math.log10(pt.getX()), pt.getY());
//				funcs.add(logHBFunc);
//				chars.add(hbChar);
//				xyzSpec.setXYElems(funcs);
//				xyzSpec.setXYChars(chars);
//			}
			
			XYZGraphPanel xyzGP = buildXYZGraphPanel();
			xyzGP.drawPlot(xyzSpec, false, false, new Range(0d, maxX+0.5*gridSpacingX),
					new Range(0d, yRange.getUpperBound()+0.5*gridSpacingY));
			// write plot
			xyzGP.getChartPanel().setSize(plotWidth, plotHeight);
			xyzGP.saveAsPNG(new File(outputDir, prefix+"_hist2D.png").getAbsolutePath());
			xyzGP.saveAsPDF(new File(outputDir, prefix+"_hist2D.pdf").getAbsolutePath());
			
//			// now write CSV
//			if (!meanSlip) {
//				IncrementalMagFreqDist magFunc = MFDPlot.buildIncrementalFunc(scatter.getMinY(), csv_mag_delta);
//				CSVFile<String> csv = new CSVFile<>(true);
//				List<String> header = new ArrayList<String>();
//				header.add("Magnitude");
//				header.add("Mean");
//				header.add("Standard Deviation");
//				for (double f : csv_fractiles)
//					header.add((float)f+" fractile");
//				csv.addLine(header);
//				
//				for (int i=0; i<magFunc.size(); i++) {
//					double mag = magFunc.getX(i);
//					double minMag = mag - 0.5*csv_mag_delta;
//					double maxMag = mag + 0.5*csv_mag_delta;
//					
//					List<Double> areasForBin = new ArrayList<>();
//					
//					for (Point2D pt : scatter)
//						if (pt.getY() >= minMag && pt.getY() < maxMag)
//							areasForBin.add(pt.getX());
//					
//					double[] areasArray = Doubles.toArray(areasForBin);
//					List<String> line = new ArrayList<>();
//					line.add((float)mag+"");
//					line.add(StatUtils.mean(areasArray)+"");
//					line.add(Math.sqrt(StatUtils.variance(areasArray))+"");
//					for (double f : csv_fractiles)
//						if (f == 0d)
//							line.add(StatUtils.min(areasArray)+"");
//						else
//							line.add(StatUtils.percentile(areasArray, f*100d)+"");
//					csv.addLine(line);
//				}
//				
//				csv.writeToFile(new File(outputDir, prefix+".csv"));
//			}
		}
		plotExample();
	}
	
	private void plotExample() throws IOException {
		if (exampleEvent == null) {
			System.out.println("No suitable example events found, skipping");
			return;
		}
		System.out.println("Plotting example rupture");
		
		List<SubSectionMapping> mappings = mapper.getAllSubSectionMappings(exampleEvent).get(0);
		
		List<XYAnnotation> anns = new ArrayList<>();
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		ArrayList<SimulatorElement> elems = exampleEvent.getAllElements();
		double[] elemSlips = exampleEvent.getAllElementSlips();
		Map<SimulatorElement, Double> slipsMap = new HashMap<>();
		for (int i=0; i<elems.size(); i++)
			slipsMap.put(elems.get(i), elemSlips[i]);
		
		CPT slipCPT = GMT_CPT_Files.GMT_HOT.instance().reverse().rescale(0d, StatUtils.max(elemSlips));
		
		double curDAS = 0;
		Map<SlipAlongSectAlgorithm, DAS_Record> algDAS_extents = new HashMap<>();
		Stroke thickElemStroke = new BasicStroke(1.5f);
		Color slippedElemColor = Color.BLACK;
		Stroke regElemStroke = new BasicStroke(1f);
		Color regElemColor = Color.GRAY;
		double maxDepth = 0d;
		List<XYAnnotation> slipAnns = new ArrayList<>();
		for (int m=0; m<mappings.size(); m++) {
			SubSectionMapping mapping = mappings.get(m);
			double totLen = mapping.getLengthForSlip(SlipAlongSectAlgorithm.MID_SEIS_FULL_LEN);
			HashSet<SimulatorElement> midSeisElems = mapper.getSlipSectionsForSect(mapping.getSubSect());
			for (SimulatorElement elem : mapper.getElementsForSection(mapping.getSubSect())) {
				SubSectDAS_Record elemDAS = mapper.getElemSubSectDAS(elem);
				DefaultXY_DataSet elemXY = new DefaultXY_DataSet();
				Vertex[] verts = elem.getVertices();
				for (int i=0; i<verts.length; i++) {
					double das = curDAS + elemDAS.vertDASs[i];
					elemXY.set(das, verts[i].getDepth());
				}
				
				boolean slipped = slipsMap.containsKey(elem);
				boolean midSeis = midSeisElems.contains(elem);
				
				Color fillColor;
				Color paint;
				if (slipped) {
					fillColor = slipCPT.getColor(slipsMap.get(elem).floatValue());
					paint = slippedElemColor;
				} else {
					fillColor = Color.WHITE;
					paint = regElemColor;
				}
				Stroke stroke;
				if (midSeis)
					stroke = thickElemStroke;
				else
					stroke = regElemStroke;
				
				if (midSeis) {
					// make it a little transparent
					fillColor = new Color(fillColor.getRed(), fillColor.getGreen(), fillColor.getBlue(), 160);
				} else {
					// make it very transparent
					fillColor = new Color(fillColor.getRed(), fillColor.getGreen(), fillColor.getBlue(), 80);
					if (!slipped)
						paint = new Color(paint.getRed(), paint.getGreen(), paint.getBlue(), 127);
				}
				
				double[] polyElems = new double[verts.length*2];
				int ind = 0;
				for (Point2D pt : elemXY) {
					polyElems[ind++] = pt.getX();
					polyElems[ind++] = pt.getY();
					maxDepth = Math.max(maxDepth, pt.getY());
				}
				XYPolygonAnnotation ann = new XYPolygonAnnotation(polyElems, stroke, paint, fillColor);
				if (slipped)
					slipAnns.add(ann);
				else
					anns.add(ann);
			}
			for (SlipAlongSectAlgorithm alg : SlipAlongSectAlgorithm.values()) {
				DAS_Record algDAS = mapping.getDASforSlip(alg);
				if (algDAS != null) {
					double startDAS = curDAS + algDAS.startDAS;
					double endDAS = curDAS + algDAS.endDAS;
					if (algDAS_extents.containsKey(alg)) {
						DAS_Record prev = algDAS_extents.get(alg);
						startDAS = Math.min(startDAS, prev.startDAS);
						endDAS = Math.max(endDAS, prev.endDAS);
					}
					algDAS_extents.put(alg, new DAS_Record(startDAS, endDAS));
				}
			}
			// draw the subsection now
			double sectUpperDepth = mapping.getSubSect().getOrigAveUpperDepth();
			double sectLowerDepth = mapping.getSubSect().getAveLowerDepth();
			maxDepth = Math.max(maxDepth, sectLowerDepth);
			DefaultXY_DataSet sectOutline = new DefaultXY_DataSet();
			sectOutline.set(curDAS, sectUpperDepth);
			sectOutline.set(curDAS+totLen, sectUpperDepth);
			sectOutline.set(curDAS+totLen, sectLowerDepth);
			sectOutline.set(curDAS, sectLowerDepth);
			sectOutline.set(curDAS, sectUpperDepth);
			if (m == 0)
				sectOutline.setName("Subsections");
			funcs.add(sectOutline);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, new Color(210, 105, 30))); // brown
			curDAS += totLen;
		}
		
		// add all slip anns on top
		anns.addAll(slipAnns);
		
		Range xRange = new Range(-curDAS*0.005, curDAS*1.1);
		
		// now draw mid-seismogenic depth range
		double[] midDepthRange = mapper.getSlipOnSectionDepthConstraints(mappings.get(0).getSubSect());
		DefaultXY_DataSet midSeisZone = new DefaultXY_DataSet();
		midSeisZone.set(xRange.getLowerBound(), midDepthRange[0]);
		midSeisZone.set(xRange.getUpperBound(), midDepthRange[0]);
		midSeisZone.set(xRange.getUpperBound(), midDepthRange[1]);
		midSeisZone.set(xRange.getLowerBound(), midDepthRange[1]);
		midSeisZone.set(xRange.getLowerBound(), midDepthRange[0]);
		midSeisZone.setName("Mid-Seismogenic Zone");
		funcs.add(midSeisZone);
		Color midColor = Color.CYAN.darker();
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, midColor));
		double midAnnX = 0.5*(xRange.getUpperBound()+curDAS);
		double midDepthSpan = midDepthRange[1]-midDepthRange[0];
		double midAnnY1 = midDepthRange[0] + 0.47*midDepthSpan;
		double midAnnY2 = midDepthRange[0] + 0.53*midDepthSpan;
		Font midFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
		XYTextAnnotation midAnn1 = new XYTextAnnotation("Mid-Seis", midAnnX, midAnnY1);
		midAnn1.setFont(midFont);
		midAnn1.setTextAnchor(TextAnchor.BASELINE_CENTER);
		midAnn1.setPaint(midColor);
		anns.add(midAnn1);
		XYTextAnnotation midAnn2 = new XYTextAnnotation("Slip Zone", midAnnX, midAnnY2);
		midAnn2.setFont(midFont);
		midAnn2.setTextAnchor(TextAnchor.TOP_CENTER);
		midAnn2.setPaint(midColor);
		anns.add(midAnn2);
		
		Color[] algColors = {Color.BLACK, Color.RED.darker(), Color.GREEN.darker(), Color.BLUE.darker(), Color.RED.darker()};
		SlipAlongSectAlgorithm[] algs = SlipAlongSectAlgorithm.values();
		double depthDelta = maxDepth*0.07;
		maxDepth += depthDelta;
		Font font = new Font(Font.SANS_SERIF, Font.BOLD, 20);
		Font subFont = new Font(Font.SANS_SERIF, Font.BOLD, 16);
		for (int i=0; i<algs.length; i++) {
			Color color = algColors[i % algColors.length];
			SlipAlongSectAlgorithm alg = algs[i];
			DAS_Record das = algDAS_extents.get(alg);
			if (das == null)
				continue;
			maxDepth += depthDelta;
			
			Stroke dotStroke = PlotLineType.DOTTED.buildStroke(3f);
			anns.add(new XYLineAnnotation(das.startDAS, midDepthRange[0], das.startDAS, maxDepth, dotStroke, color));
			anns.add(new XYLineAnnotation(das.endDAS, midDepthRange[0], das.endDAS, maxDepth, dotStroke, color));

			Stroke regStroke = new BasicStroke(3f);
			anns.add(new XYLineAnnotation(das.startDAS, maxDepth, das.endDAS, maxDepth, regStroke, color));
			Stroke thickStroke = new BasicStroke(5f);
			anns.add(new XYLineAnnotation(das.startDAS, midDepthRange[0], das.startDAS, midDepthRange[1], thickStroke, color));
			anns.add(new XYLineAnnotation(das.endDAS, midDepthRange[0], das.endDAS, midDepthRange[1], thickStroke, color));

			double sumLength = 0d;
			double areaWeightedSlip = 0d;
			double sumArea = 0d;
			for (SubSectionMapping mapping : mappings) {
				sumLength += mapping.getLengthForSlip(alg);
				double area = mapping.getAreaForAverageSlip(alg);
				double slip = mapping.getAverageSlip(alg);
				areaWeightedSlip += slip*area;
				sumArea += area;
			}
			double aveSlip = areaWeightedSlip == 0 ? 0 : areaWeightedSlip / sumArea;
			
			XYTextAnnotation textAnn = new XYTextAnnotation(alg+"", das.midDAS, maxDepth-0.7*depthDelta);
			textAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
			textAnn.setFont(font);
			textAnn.setPaint(color);
			anns.add(textAnn);
			String slipStr = "len="+optionalDigitDF.format(sumLength)+" [km], ⟨slip⟩="+optionalDigitDF.format(aveSlip)+" [m]";
			textAnn = new XYTextAnnotation(slipStr, das.midDAS, maxDepth-0.1*depthDelta);
			textAnn.setTextAnchor(TextAnchor.BASELINE_CENTER);
			textAnn.setFont(subFont);
			textAnn.setPaint(color);
			anns.add(textAnn);
			maxDepth += depthDelta;
		}
		
		Range yRange = new Range(maxDepth*-0.005, maxDepth);
		
		String title = "Slip-Length Algorithms Example";
		String xAxisLabel = "Distance Along Strike (km)";
		String yAxisLabel = "Depth (km)";
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		plot.setLegendVisible(false);
		plot.setPlotAnnotations(anns);
		
		double xPerY = xRange.getLength()/yRange.getLength();
		int ySize = 800;
		int xSize = (int)Math.round(ySize*xPerY);
		
		HeadlessGraphPanel gp = buildGraphPanel();
		gp.setyAxisInverted(true);
		gp.drawGraphPanel(plot, false, false, xRange, yRange);
		gp.getChartPanel().setSize(xSize, ySize);
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix()+"_example_rupture";
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		
		double skipYears = 5000;
		double minMag = 6.5;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		SlipLengthScalingPlot plot = new SlipLengthScalingPlot(catalog.getSubSectMapper(), minMag);
		plot.initialize(catalog.getName(), outputDir, "slip_len");
		
		EXAMPLE_DEBUG = true;
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
	}

}
