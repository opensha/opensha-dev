package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.DAS_Record;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SlipAlongSectAlgorithm;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SlipAlongRupturePlot extends AbstractPlot {
	
	private static final boolean D = true;
	private static final boolean DD = false;
	
	private RSQSimSubSectionMapper mapper;
	private double minMag;
	private SlipAlongSectAlgorithm slipAlg;
	
	private Map<Integer, HashSet<Integer>> faultParentIDsMap;
	
	private EvenlyDiscrXYZ_DataSet[] normalizedSingleFaultXYZ;
	private EvenlyDiscrXYZ_DataSet[] absoluteSingleFaultXYZ;
	private EvenlyDiscrXYZ_DataSet[] normalizedMultiFaultXYZ;
	private EvenlyDiscrXYZ_DataSet[] absoluteMultiFaultXYZ;
	private List<Range> lengthBins;

	public SlipAlongRupturePlot(RSQSimSubSectionMapper mapper, double minMag, SlipAlongSectAlgorithm slipAlg, FaultModels fm,
			List<Range> lengthBins) {
		this.mapper = mapper;
		this.minMag = minMag;
		this.slipAlg = slipAlg;
		this.lengthBins = lengthBins;
		
		Map<String, List<Integer>> fmNamesMap = fm.getNamedFaultsMapAlt();
		faultParentIDsMap = new HashMap<>();
		for (List<Integer> parentIDs : fmNamesMap.values()) {
			HashSet<Integer> set = new HashSet<>(parentIDs);
			for (Integer parent : parentIDs)
				faultParentIDsMap.put(parent, set);
		}
		
		mapper.trackSlipOnSections();
		
		int numY = 20;
		double minY = 0.025;
		double deltaY = 0.05;
		
		int numX = 20;
		double deltaNormX = 0.5d/(double)numX;
		double minNormX = deltaNormX*0.5;
		
		normalizedSingleFaultXYZ = new EvenlyDiscrXYZ_DataSet[lengthBins.size()];
		absoluteSingleFaultXYZ = new EvenlyDiscrXYZ_DataSet[lengthBins.size()];
		normalizedMultiFaultXYZ = new EvenlyDiscrXYZ_DataSet[lengthBins.size()];
		absoluteMultiFaultXYZ = new EvenlyDiscrXYZ_DataSet[lengthBins.size()];
		
		for (int i=0; i<lengthBins.size(); i++) {
			Range lengthBin = lengthBins.get(i);
			double maxAbsX = lengthBin == null || Double.isInfinite(lengthBin.getUpperBound()) ? 200 : lengthBin.getUpperBound()/2;
			double deltaAbsX = maxAbsX/(double)numX;
			double minAbsX = deltaAbsX*0.5;
			normalizedSingleFaultXYZ[i] = new EvenlyDiscrXYZ_DataSet(numX, numY, minNormX, minY, deltaNormX, deltaY);
			absoluteSingleFaultXYZ[i] = new EvenlyDiscrXYZ_DataSet(numX, numY, minAbsX, minY, deltaAbsX, deltaY);
			normalizedMultiFaultXYZ[i] = new EvenlyDiscrXYZ_DataSet(numX*2, numY, -normalizedSingleFaultXYZ[i].getMaxX(), minY, deltaNormX, deltaY);
			absoluteMultiFaultXYZ[i] = new EvenlyDiscrXYZ_DataSet(numX*2, numY, -absoluteSingleFaultXYZ[i].getMaxX(), minY, deltaAbsX, deltaY);
		}
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		if (e.getMagnitude() < minMag)
			return;
		Preconditions.checkState(e instanceof RSQSimEvent);
		List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings((RSQSimEvent)e);
		
		// now bundle by fault
		List<List<SubSectionMapping>> faultBundledMappings = new ArrayList<>();
		Integer curParentID = null;
		List<SubSectionMapping> curFaultMappings = null;
		double maxSlip = 0d;
		List<Double> indvMaxSlips = new ArrayList<>();
		for (List<SubSectionMapping> parentMappings : mappings) {
			int parentID = parentMappings.get(0).getSubSect().getParentSectionId();
			if (curParentID == null) {
				// first one
				curParentID = parentID;
				curFaultMappings = new ArrayList<>();
				faultBundledMappings.add(curFaultMappings);
				indvMaxSlips.add(0d);
			} else {
				// let's see if this parent is also part of that fault
				boolean match = faultParentIDsMap.containsKey(parentID) && faultParentIDsMap.get(parentID).contains(curParentID);
				if (!match) {
					// it's a new fault
					curParentID = parentID;
					curFaultMappings = new ArrayList<>();
					faultBundledMappings.add(curFaultMappings);
					indvMaxSlips.add(0d);
				}
			}
			curFaultMappings.addAll(parentMappings);
			double myMaxSlip = 0d;
			for (SubSectionMapping mapping : parentMappings)
				myMaxSlip = Math.max(myMaxSlip, mapping.getAverageSlip(slipAlg));
			maxSlip = Math.max(maxSlip, myMaxSlip);
			indvMaxSlips.set(indvMaxSlips.size()-1, Math.max(indvMaxSlips.get(indvMaxSlips.size()-1), myMaxSlip));
		}
		
		if (faultBundledMappings.size() == 1) {
			List<SubSectionMapping> bundle = faultBundledMappings.get(0);
			if (bundle.size() > 1)
				processFaultBundle(bundle, false, false, maxSlip, absoluteSingleFaultXYZ, normalizedSingleFaultXYZ);
		} else {
			SubSectionMapping firstMapping = mappings.get(0).get(0);
			SubSectionMapping lastMapping = mappings.get(mappings.size()-1).get(mappings.get(mappings.size()-1).size()-1);
			FaultTrace firstTrace = firstMapping.getSubSect().getFaultTrace();
			FaultTrace lastTrace = lastMapping.getSubSect().getFaultTrace();
			Location firstLoc = firstMapping.isReversed() ? firstTrace.last() : firstTrace.first();
			Location lastLoc = lastMapping.isReversed() ? lastTrace.last() : lastTrace.first();
			for (int i=1; i<faultBundledMappings.size(); i++) {
				List<SubSectionMapping> bundle0 = faultBundledMappings.get(i-1);
				List<SubSectionMapping> bundle1 = faultBundledMappings.get(i);
				if (bundle0.size() > 1 && bundle1.size() > 1) {
					// figure out which one should be on the left vs right
					// first determine the midpoint of the jump
					SubSectionMapping beforeMapping = bundle0.get(bundle0.size() - 1);
					SubSectionMapping afterMapping = bundle1.get(0);
					FaultTrace beforeJumpTrace = beforeMapping.getSubSect().getFaultTrace();
					FaultTrace afterJumpTrace = afterMapping.getSubSect().getFaultTrace();
					Location beforeJumpLoc = beforeMapping.isReversed() ? beforeJumpTrace.last() : beforeJumpTrace.first();
					Location afterJumpLoc = afterMapping.isReversed() ? afterJumpTrace.last() : afterJumpTrace.first();
					LocationVector vector = LocationUtils.vector(beforeJumpLoc, afterJumpLoc);
					vector.setHorzDistance(vector.getHorzDistance()*0.5);
					Location midJumpLoc = LocationUtils.location(beforeJumpLoc, vector);

					double distToStart = LocationUtils.horzDistanceFast(midJumpLoc, firstLoc);
					double distToEnd = LocationUtils.horzDistanceFast(midJumpLoc, lastLoc);
					boolean beforeOnLeft = distToStart < distToEnd;
					
					double slip0 = indvMaxSlips.get(i-1);
					double slip1 = indvMaxSlips.get(i);
					double myMaxSlip = Math.max(slip0, slip1);
					processFaultBundle(bundle0, true, beforeOnLeft, myMaxSlip,
							absoluteMultiFaultXYZ, normalizedMultiFaultXYZ);
					processFaultBundle(bundle1, true, !beforeOnLeft, myMaxSlip,
							absoluteMultiFaultXYZ, normalizedMultiFaultXYZ);
				}
			}
		}
	}
	
	private void processFaultBundle(List<SubSectionMapping> faultMappings, boolean multiFault, boolean leftSide, double maxSlip,
			EvenlyDiscrXYZ_DataSet[] absXYZ, EvenlyDiscrXYZ_DataSet[] normXYZ) {
		double curX = 0;
		DefaultXY_DataSet myFunc = new DefaultXY_DataSet();
		for (int i=0; i<faultMappings.size(); i++) {
			SubSectionMapping mapping = faultMappings.get(i);
			DAS_Record slipDAS = mapping.getDASforSlip(slipAlg);
			
			FaultTrace trace = mapping.getSubSect().getFaultTrace();
			double len = trace.getTraceLength();
			
			if (slipDAS != null) {
				double aveSlip = mapping.getAverageSlip(slipAlg);
				double normSlip = aveSlip / maxSlip;
				if (slipDAS.endDAS > len) {
					double diff = slipDAS.endDAS - len;
					Preconditions.checkState(diff < 0.01, "endDAS=% is greater than fault length=%s", slipDAS.endDAS, len);
					slipDAS = new DAS_Record(slipDAS.startDAS, len);
				}
				if (mapping.isReversed())
					slipDAS = slipDAS.getReversed(len);
				Preconditions.checkState(slipDAS.startDAS >= 0f);
				Preconditions.checkState(slipDAS.endDAS >= 0f);
				if (aveSlip > 0 && slipDAS.endDAS > slipDAS.startDAS) {
					double start, end;
					if (i == 0) {
						// first one, start where it starts not start of fault
						start = 0d;
						end = slipDAS.endDAS - slipDAS.startDAS;
						curX = len - slipDAS.startDAS;
					} else if (i == faultMappings.size()-1) {
						// last one, end where it ends
						start = curX + slipDAS.startDAS;
						end = curX + slipDAS.endDAS;
						curX += slipDAS.endDAS;
					} else {
						start = curX + slipDAS.startDAS;
						end = curX + slipDAS.endDAS;
						curX += len;
					}
	
					myFunc.set(start, normSlip);
					myFunc.set(end, normSlip);
				}
			}
			
			if (DD) {
				System.out.println("i="+i+", curX="+curX+", maxX="+myFunc.getMaxX()+", reversed="+mapping.isReversed());
				if (slipDAS != null)
					System.out.println("\tslipDAS=["+slipDAS.startDAS+" "+slipDAS.endDAS+"]");
			}
		}
		double totalLen = curX;
		if (myFunc.size() == 0)
			return;
		
		validateFuncOrder(myFunc);

		final boolean debugPlot = DD && Math.random() < 0.0001 && numDebugPlots < 10 && multiFault;
		for (boolean normalized : new boolean[] {false, true}) {
			XY_DataSet func = normalized ? getNormalized(myFunc) : myFunc;
			XY_DataSet firstHalf = getFirstHalf(func);
			if (leftSide)
				firstHalf = getMirroredNegative(firstHalf);
			
			XY_DataSet mirrored = null; 
			if (!multiFault)
				// use both ends of this fault
				mirrored = getFirstHalf(getMirroredPositive(func));
			
			EvenlyDiscrXYZ_DataSet[] xyz = normalized ? normXYZ : absXYZ;
			for (int i=0; i<lengthBins.size(); i++) {
				Range lengthBin = lengthBins.get(i);
				if (lengthBin != null && !lengthBin.contains(totalLen))
					continue;
				processSlipFunc(firstHalf, xyz[i]);
				if (mirrored != null)
					processSlipFunc(mirrored, xyz[i]);
				if (debugPlot && lengthBin != null)
					debugPlot(xyz[i], firstHalf, mirrored);
			}
		}
	}
	
	private static int numDebugPlots = 0;
	private void debugPlot(EvenlyDiscrXYZ_DataSet refXYZ, XY_DataSet forward, XY_DataSet mirrored) {
		// create an empty copy
		EvenlyDiscrXYZ_DataSet xyz = refXYZ.copy();
		for (int i=0; i<xyz.size(); i++)
			xyz.set(i, 0d);
		
		List<XY_DataSet> xyFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> xyChars = new ArrayList<>();
		
		processSlipFunc(forward, xyz);
		xyFuncs.add(forward);
		xyChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		if (mirrored != null) {
			processSlipFunc(mirrored, xyz);
			xyFuncs.add(mirrored);
			xyChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		}
		
		double minNonZero = Double.POSITIVE_INFINITY;
		for (int i=0; i<xyz.size(); i++) {
			double val = xyz.get(i);
			if (val > 0 && val < minNonZero)
				minNonZero = val;
		}
		
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		if (!Double.isInfinite(minNonZero)) {
			Range logRange;
			if (minNonZero == xyz.getMaxZ())
				logRange = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))), Math.pow(10, Math.ceil(Math.log10(minNonZero)+0.001)));
			else
				logRange = calcEncompassingLog10Range(minNonZero, xyz.getMaxZ());
			cpt = cpt.rescale(Math.log10(logRange.getLowerBound()), Math.log10(logRange.getUpperBound()));
		}
		xyz.log10();
		
		XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, "Debug Event", "X", "Y", null);
		xyzSpec.setXYElems(xyFuncs);
		xyzSpec.setXYChars(xyChars);

		double minX = xyz.getMinX();
		double maxX = xyz.getMaxX();
		double gridSpacingX = xyz.getGridSpacingX();
		double minY = xyz.getMinY();
		double maxY = xyz.getMaxY();
		double gridSpacingY = xyz.getGridSpacingY();
		
		XYZGraphPanel xyzGP = buildXYZGraphPanel();
		xyzGP.drawPlot(xyzSpec, false, false, new Range(minX-0.5*gridSpacingX, maxX+0.5*gridSpacingX),
				new Range(minY-0.5*gridSpacingY, maxY+0.5*gridSpacingY));
		
		new XYZPlotWindow(xyzGP);
		numDebugPlots++;
	}
	
	private static void validateFuncOrder(XY_DataSet func) {
		for (int i=1; i<func.size(); i++) {
			double prevX = func.getX(i-1);
			double x = func.getX(i);
			Preconditions.checkState(x >= prevX, "Dataset out of order? x[%s]=%s, x[%s]=%s\n%s", i-1, prevX, i, x, func);
		}
	}
	
	private static XY_DataSet getNormalized(XY_DataSet slipFunc) {
		if (D) validateFuncOrder(slipFunc);
		double minX = slipFunc.getMinX();
		double maxX = slipFunc.getMaxX();
		double spanX = maxX - minX;
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (int i=0; i<slipFunc.size(); i+=2) {
			double start = (slipFunc.getX(i)-minX)/spanX;
			double end = (slipFunc.getX(i+1)-minX)/spanX;
			double normSlip = slipFunc.getY(i);
			ret.set(start, normSlip);
			ret.set(end, normSlip);
		}
		return ret;
	}
	
	private static XY_DataSet getFirstHalf(XY_DataSet slipFunc) {
		if (D) validateFuncOrder(slipFunc);
		double minX = slipFunc.getMinX();
		double maxX = slipFunc.getMaxX();
		Preconditions.checkState(slipFunc.size() >= 2);
		Preconditions.checkState(Double.isFinite(minX));
		Preconditions.checkState(minX >= 0d);
		Preconditions.checkState(maxX >= minX);
		double midX = 0.5*(minX + maxX);
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (int i=0; i<slipFunc.size(); i+=2) {
			double start = slipFunc.getX(i);
			if (start > midX)
				break;
			double end = Math.min(midX, slipFunc.getX(i+1));
			double normSlip = slipFunc.getY(i);
			ret.set(start, normSlip);
			ret.set(end, normSlip);
		}
		Preconditions.checkState(ret.size() >= 2, "Error getting first half. MinX=%s, MaxX=%s, MidX=%s, size=%s, first=%s",
				minX, maxX, midX, ret.size(), slipFunc.getX(0));
		return ret;
	}
	
	private static XY_DataSet getMirroredPositive(XY_DataSet slipFunc) {
		if (D) validateFuncOrder(slipFunc);
		double maxX = slipFunc.getMaxX();
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (int i=slipFunc.size()-1; i>0; i-=2) {
			double start = slipFunc.getX(i-1);
			double end = slipFunc.getX(i);
			double normSlip = slipFunc.getY(i);
			ret.set(maxX - end, normSlip);
			ret.set(maxX - start, normSlip);
		}
		return ret;
	}
	
	private static XY_DataSet getMirroredNegative(XY_DataSet slipFunc) {
		if (D) validateFuncOrder(slipFunc);
		double minX = slipFunc.getMinX();
		Preconditions.checkState(minX >= 0d, "Mirroring to be negative, but already have a negative minX=%s", minX);
		DefaultXY_DataSet ret = new DefaultXY_DataSet();
		for (int i=slipFunc.size()-1; i>0; i-=2) {
			double start = slipFunc.getX(i-1);
			double end = slipFunc.getX(i);
			double normSlip = slipFunc.getY(i);
			ret.set(-end, normSlip);
			ret.set(-start, normSlip);
		}
		return ret;
	}
	
	private static void processSlipFunc(XY_DataSet slipFunc, EvenlyDiscrXYZ_DataSet xyz) {
		double binWidth = xyz.getGridSpacingX();
		double halfBin = binWidth*0.5;
		if (D) validateFuncOrder(slipFunc);
		for (int i=0; i<slipFunc.size(); i+=2) {
			double start = slipFunc.getX(i);
			double end = slipFunc.getX(i+1);
			double normSlip = slipFunc.getY(i);
			int xInd = xyz.getXIndex(start);
			int yInd = xyz.getYIndex(normSlip);
			while (xInd < xyz.getNumX()) {
				if (xInd < 0) {
					xInd++;
					continue;
				}
				double binMiddle = xyz.getX(xInd);
				double binStart = binMiddle - halfBin;
				double binEnd = binMiddle + halfBin;
				if (end < binStart)
					break;
				if (binEnd < start) {
					xInd++;
					continue;
				}
				double overlapStart = Math.max(start, binStart);
				double overlapEnd = Math.min(end, binEnd);
				Preconditions.checkState(overlapEnd >= overlapStart);
				double overlapSpan = overlapEnd - overlapStart;
				Preconditions.checkState((float)overlapSpan <= (float)binWidth, "overlap span is %s (%s - %s) but width is %s for bin [%s %s]",
						overlapSpan, overlapEnd, overlapStart, binWidth, binStart, binEnd);
				double fract = overlapSpan / binWidth;
				
				xyz.set(xInd, yInd, xyz.get(xInd, yInd) + fract);
				xInd++;
			}
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		String name = getCatalogName();
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix();

		for (int i=0; i<lengthBins.size(); i++) {
			Range lengthBin = lengthBins.get(i);
			String lenStr;
			String prefixAdd = "";
			if (lengthBin == null) {
				lenStr = ", All Lengths";
			} else {
				if (Double.isInfinite(lengthBin.getUpperBound())) {
					lenStr = ", Len≥"+optionalDigitDF.format(lengthBin.getLowerBound());
					prefixAdd = "_len_"+optionalDigitDF.format(lengthBin.getLowerBound())+"+";
				} else {
					lenStr = ", Len=["+optionalDigitDF.format(lengthBin.getLowerBound())+" "+optionalDigitDF.format(lengthBin.getUpperBound())+"]";
					prefixAdd = "_len_"+optionalDigitDF.format(lengthBin.getLowerBound())+"_"+optionalDigitDF.format(lengthBin.getUpperBound());
				}
			}
			writePlot(outputDir, prefix+"_single_abs"+prefixAdd, "Single Fault Slip Along Rupture"+lenStr, "Distance Along Strike (km)",
					absoluteSingleFaultXYZ[i]);
			writePlot(outputDir, prefix+"_single_norm"+prefixAdd, "Single Fault Slip Along Rupture"+lenStr, "Normalized Distance Along Strike",
					normalizedSingleFaultXYZ[i]);
			writePlot(outputDir, prefix+"_multi_abs"+prefixAdd, "Multi-Fault Slip Along Rupture"+lenStr, "Distance Along Strike (km)",
					absoluteMultiFaultXYZ[i]);
			writePlot(outputDir, prefix+"_multi_norm"+prefixAdd, "Multi-Fault Slip Along Rupture"+lenStr, "Normalized Distance Along Strike",
					normalizedMultiFaultXYZ[i]);
		}
	}
	
	private void writePlot(File outputDir, String prefix, String title, String xAxisLabel, EvenlyDiscrXYZ_DataSet xyz) throws IOException {
		String yAxisLabel = "Normalized Slip";
		String zAxisLabel = "Log10(Density)";
		
		ArbitrarilyDiscretizedFunc meanRightFunc = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc meanLeftFunc = xyz.getMinX() < 0 ? new ArbitrarilyDiscretizedFunc() : null;
		ArbitrarilyDiscretizedFunc sqrtSinFunc = xyz.getMaxX() > 1 ? null : new ArbitrarilyDiscretizedFunc();
		
		double maxLeft = 0d;
		double maxRight = 0d;
		
//		if (sqrtSinFunc != null)
//			System.out.println("Doing "+prefix);
		for (int x=0; x<xyz.getNumX(); x++) {
			double xVal = xyz.getX(x);
			boolean left = xVal < 0;
			double meanY = 0d;
			double sumWeights = 0d;
			for (int y=0; y<xyz.getNumY(); y++) {
				double weight = xyz.get(x, y);
				sumWeights += weight;
				meanY += xyz.getY(y)*weight;
				if (left)
					maxLeft = Math.max(weight, maxLeft);
				else
					maxRight = Math.max(weight, maxRight);
			}
			meanY /= sumWeights;
			if (meanY > 0) {
				if (xVal < 0)
					meanLeftFunc.set(xVal, meanY);
				else
					meanRightFunc.set(xVal, meanY);
			}
			if (sqrtSinFunc != null) {
//				System.out.println("\tsumY["+(float)xVal+"]: "+sumWeights);
				double relX = xVal * Math.PI;
				sqrtSinFunc.set(xVal, Math.sqrt(Math.abs(Math.sin(relX))));
			}
		}
		if (sqrtSinFunc != null)
			// fill in zero val
			sqrtSinFunc.set(0d, 0d);
		
		xyz = xyz.copy();
		
		// now normalize. if this is a multi-fault plot, normalize each side separately
		// (can have different amounts on either side)
		
		for (int xInd=0; xInd<xyz.getNumX(); xInd++) {
			double scale = xyz.getX(xInd) < 0 ? maxLeft : maxRight;
			for (int yInd=0; yInd<xyz.getNumY(); yInd++)
				xyz.set(xInd, yInd, xyz.get(xInd, yInd)/scale);
		}
		
		double minNonZero = Double.POSITIVE_INFINITY;
		for (int i=0; i<xyz.size(); i++) {
			double val = xyz.get(i);
			if (val > 0 && val < minNonZero)
				minNonZero = val;
		}
		Range logRange = calcEncompassingLog10Range(Math.max(minNonZero, 1d/Math.max(maxLeft, maxRight)), xyz.getMaxZ());
		xyz.log10();
		for (int i=0; i<xyz.size(); i++)
			if (Double.isInfinite(xyz.get(i)))
				xyz.set(i, Double.NaN);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(Math.log10(logRange.getLowerBound()), Math.log10(logRange.getUpperBound()));
		cpt.setNanColor(new Color(255, 255, 255, 0));
		
		XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, title, xAxisLabel, yAxisLabel, zAxisLabel);
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		meanRightFunc.setName("Mean");
		funcs.add(meanRightFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		if (meanLeftFunc != null) {
			funcs.add(meanLeftFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		}
		
		if (sqrtSinFunc != null) {
			sqrtSinFunc.setName("sqrt(|sin(x*π)|)");
			funcs.add(sqrtSinFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		}
		
		xyzSpec.setXYElems(funcs);
		xyzSpec.setXYChars(chars);

		double minX = xyz.getMinX();
		double maxX = xyz.getMaxX();
		double gridSpacingX = xyz.getGridSpacingX();
		double minY = xyz.getMinY();
		double maxY = xyz.getMaxY();
		double gridSpacingY = xyz.getGridSpacingY();
		
		XYZGraphPanel xyzGP = buildXYZGraphPanel();
		xyzGP.drawPlot(xyzSpec, false, false, new Range(minX-0.5*gridSpacingX, maxX+0.5*gridSpacingX),
				new Range(minY-0.5*gridSpacingY, maxY+0.5*gridSpacingY));
		// write plot
		xyzGP.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
		xyzGP.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		xyzGP.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}
	
	public static void main(String[] args) {
		try {
			File baseDir = new File("/home/kevin/Simulators/catalogs");
			
			double skipYears = 5000;
			double minMag = 6.5;
			
			RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
			
			File outputDir = new File("/tmp");
			
			List<Range> lengthBins = new ArrayList<>();
			lengthBins.add(null);
			lengthBins.add(new Range(0d, 20));
			lengthBins.add(new Range(20, 50));
			lengthBins.add(new Range(50, 100));
			lengthBins.add(new Range(100, Double.POSITIVE_INFINITY));
			
			SlipAlongRupturePlot plot = new SlipAlongRupturePlot(catalog.getSubSectMapper(), minMag,
					SlipAlongSectAlgorithm.MID_SEIS_SLIPPED_LEN, FaultModels.FM3_1, lengthBins);
			plot.initialize(catalog.getName(), outputDir, "slip_along_rup");
			
			for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
				plot.processEvent(e);
			
			System.out.println("Finalizing plot...");
			
			plot.finalizePlot();

			System.out.println("DONE");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

}
