package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.RSQSimEventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.aftershockStatistics.AftershockStatsCalc;

public class RSQSimPatchBvalCalc {

	public static void main(String[] args) throws IOException {
		File dir = new File("/data/rsqsim/supraSeisGeo2");
		File nuclOutputFile = new File("/tmp/nucl_b_vals.txt");
		File partOutputFile = new File("/tmp/part_b_vals.txt");
		File geomFile = new File(dir, "UCERF3.D3.1.1km.tri.2.flt");
		
		List<SimulatorElement> elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
		int patchOffset = elements.get(0).getID();
		System.out.println("Patch offset: "+patchOffset);
		
		double minMag = 5;
		List<RuptureIdentifier> loadIdens = Lists.newArrayList();
		loadIdens.add(new MagRangeRuptureIdentifier(minMag, Double.POSITIVE_INFINITY));
		
		IncrementalMagFreqDist[] patchNuclMFDs = new IncrementalMagFreqDist[elements.size()];
		IncrementalMagFreqDist[] patchPartMFDs = new IncrementalMagFreqDist[elements.size()];
		double mfdDelta = 0.1;
		double mfdMin = minMag + 0.5*mfdDelta;
		int mfdNum = (int)((9d - minMag)/mfdDelta);
		for (int i=0; i<elements.size(); i++) {
			patchNuclMFDs[i] = new IncrementalMagFreqDist(mfdMin, mfdNum, mfdDelta);
			patchPartMFDs[i] = new IncrementalMagFreqDist(mfdMin, mfdNum, mfdDelta);
		}
		
		Map<Integer, Double> sectMinDepthMap = Maps.newHashMap();
		Map<Integer, Double> sectMaxDepthMap = Maps.newHashMap();
		
		for (SimulatorElement e : elements) {
			Integer sectID = e.getSectionID();
			double depth = e.getAveDepth();
			
			double minDepth, maxDepth;
			if (sectMinDepthMap.containsKey(sectID)) {
				minDepth = Math.min(depth, sectMinDepthMap.get(sectID));
				maxDepth = Math.max(depth, sectMaxDepthMap.get(sectID));
			} else {
				minDepth = depth;
				maxDepth = depth;
			}
			
			sectMinDepthMap.put(sectID, minDepth);
			sectMaxDepthMap.put(sectID, maxDepth);
		}
		
		int numDepth = 10;
		IncrementalMagFreqDist[] depthMFDs = new IncrementalMagFreqDist[numDepth];
		for (int i=0; i<numDepth; i++)
			depthMFDs[i] = new IncrementalMagFreqDist(mfdMin, mfdNum, mfdDelta);
		
		System.out.println("Getting iterable");
		Iterable<RSQSimEvent> iterable = RSQSimFileReader.getEventsIterable(dir, elements, loadIdens);
		System.out.println("Got iterable");
		
		int count = 0;
		for (RSQSimEvent e : iterable) {
			if (count % 100000 == 0)
				System.out.println("working on event "+count);
			count++;
			RSQSimEventRecord rec = (RSQSimEventRecord)e.get(0);
			int firstPatch = rec.getFirstPatchToSlip();
			double mag = e.getMagnitude();
			if (mag < minMag)
				continue;
			Preconditions.checkState(firstPatch >= 0);
			int magIndex = patchNuclMFDs[firstPatch-patchOffset].getClosestXIndex(mag);
			patchNuclMFDs[firstPatch-patchOffset].add(magIndex, 1d);
			for (int patchID : e.getAllElementIDs())
				patchPartMFDs[patchID-patchOffset].add(magIndex, 1d);
			
			// now handle depth MFDs
			SimulatorElement elem = elements.get(firstPatch-patchOffset);
			double depth = elem.getAveDepth();
			Integer sectID = elem.getSectionID();
			double sectMinDepth = sectMinDepthMap.get(sectID);
			double sectMaxDepth = sectMaxDepthMap.get(sectID);
			double fractDepth = (depth - sectMinDepth)/(sectMaxDepth - sectMinDepth);
			Preconditions.checkState(fractDepth >= 0 && fractDepth <= 1d, "Bad fract depth: %s", fractDepth);
			for (int i=0; i<numDepth; i++) {
				double myVal = (double)(i+1)/(double)numDepth;
				if (fractDepth > myVal)
					continue;
				depthMFDs[i].add(magIndex, 1d);
			}
		}
		System.out.println("Processed "+count+" events");
		
		// calculate b values
		FileWriter partFW = new FileWriter(partOutputFile);
		FileWriter nuclFW = new FileWriter(nuclOutputFile);
		
		for (int i=0; i<elements.size(); i++) {
			nuclFW.write(getB(patchNuclMFDs[i])+"\n");
			partFW.write(getB(patchPartMFDs[i])+"\n");
		}
		
		partFW.close();
		nuclFW.close();
		
		// now plot depth specific
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		CPT depthCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		double maxY = 0;
		for (int i=0; i<numDepth; i++) {
			double myVal = (double)(i+1)/(double)numDepth;
			Color c = depthCPT.getColor((float)myVal);
			
			if (i == 0)
				depthMFDs[i].setName("Fract depth <= "+(float)myVal);
			else
				depthMFDs[i].setName("<= "+(float)myVal);
			maxY = Math.max(maxY, depthMFDs[i].getMaxY());
			funcs.add(depthMFDs[i]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Nucl Depth Dependent MFDs", "Magnitude", "Number");
		spec.setLegendVisible(true);
		GraphWindow gw = new GraphWindow(spec);
		gw.setXLog(false);
		gw.setYLog(true);
		gw.setAxisRange(minMag, 9d, 1d, maxY);
		
		System.out.println("DONE");
//		System.exit(0);
	}
	
	private static double getB(IncrementalMagFreqDist mfd) {
		double sumY = mfd.calcSumOfY_Vals();
		if (sumY == 0d)
			return Double.NaN;
		double magPrecision = mfd.getDelta();
		double magComplete = mfd.getMinX() - 0.5*mfd.getDelta();
		double magMean = 0;
		for (int i=0; i<mfd.size(); i++)
			magMean += mfd.getX(i)*mfd.getY(i);
		magMean /= mfd.calcSumOfY_Vals();
		return AftershockStatsCalc.getMaxLikelihood_b_value(magMean, magComplete, magPrecision);
	}

}
