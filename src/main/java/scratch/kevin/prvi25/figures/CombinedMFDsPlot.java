package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Direct;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Exact;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.latex.LaTeXUtils;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class CombinedMFDsPlot {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution crustalSol = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED);
		
		FaultSystemSolution subductionSol1 = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		FaultSystemSolution subductionSol2 = FaultSystemSolution.load(SUBDUCTION_SOL_SMALL);
		
		File outputDir = new File(FIGURES_DIR, "combined_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(9.45);
		
		List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
		List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		List<EvenlyDiscretizedFunc> texFuncs = new ArrayList<>();
		List<String> texPrefixes = new ArrayList<>();
		
		Region crustalReg = PRVI25_SeismicityRegions.CRUSTAL.load();
		
		IncrementalMagFreqDist crustalFault = calcFaultMFD(crustalReg, crustalSol, refMFD);
		crustalFault.setName("Crustal, Fault");
		incrFuncs.add(crustalFault);
		cmlFuncs.add(crustalFault.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("CrustalFault");
		
		IncrementalMagFreqDist crustalGrid = calcGriddedMFD(crustalReg, TectonicRegionType.ACTIVE_SHALLOW, crustalSol, refMFD);
		crustalGrid.setName("Crustal, Gridded");
		incrFuncs.add(crustalGrid);
		cmlFuncs.add(crustalGrid.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightblue));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("CrustalGridded");
		
		texFuncs.add(sum(crustalFault, crustalGrid).getCumRateDistWithOffset());
		texPrefixes.add("CrustalTotal");
		
		IncrementalMagFreqDist faultMuertos = average(
				calcFaultMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), subductionSol1, refMFD),
				calcFaultMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), subductionSol2, refMFD));
		faultMuertos.setName("Muertos, Interface Fault");
		incrFuncs.add(faultMuertos);
		cmlFuncs.add(faultMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.green));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("MuertosInterfaceFault");
		
		IncrementalMagFreqDist griddedMuertos = average(
				calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol1, refMFD),
				calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol2, refMFD));
		griddedMuertos.setName("Muertos, Interface Gridded");
		incrFuncs.add(griddedMuertos);
		cmlFuncs.add(griddedMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightgreen));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("MuertosInterfaceGridded");
		
		texFuncs.add(sum(faultMuertos, griddedMuertos).getCumRateDistWithOffset());
		texPrefixes.add("MuertosInterfaceTotal");
		
		IncrementalMagFreqDist slabMuertos = calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTRASLAB.load(),
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol1, refMFD);
		slabMuertos.setName("Muertos, Slab");
		incrFuncs.add(slabMuertos);
		cmlFuncs.add(slabMuertos.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightgreen));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("MuertosSlab");
		
		IncrementalMagFreqDist faultCar = average(
				calcFaultMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), subductionSol1, refMFD),
				calcFaultMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), subductionSol2, refMFD));
		faultCar.setName("Caribbean, Interface Fault");
		incrFuncs.add(faultCar);
		cmlFuncs.add(faultCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("CarInterfaceFault");
		
		IncrementalMagFreqDist griddedCar = average(
				calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol1, refMFD),
				calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTERFACE.load(), TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol2, refMFD));
		griddedCar.setName("Caribbean, Interface Gridded");
		incrFuncs.add(griddedCar);
		cmlFuncs.add(griddedCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightorange));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("CarInterfaceGridded");
		
		texFuncs.add(sum(faultCar, griddedCar).getCumRateDistWithOffset());
		texPrefixes.add("CarInterfaceTotal");
		
		IncrementalMagFreqDist slabCar = calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTRASLAB.load(),
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol1, refMFD);
		slabCar.setName("Caribbean, Slab");
		incrFuncs.add(slabCar);
		cmlFuncs.add(slabCar.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightorange));
		texFuncs.add(cmlFuncs.get(cmlFuncs.size()-1));
		texPrefixes.add("CarSlab");
		
		texFuncs.add(sum(faultCar, griddedCar, faultMuertos, griddedMuertos).getCumRateDistWithOffset());
		texPrefixes.add("InterfaceTotal");
		
		texFuncs.add(sum(slabCar, slabMuertos).getCumRateDistWithOffset());
		texPrefixes.add("SlabTotal");
		
		SummedMagFreqDist sum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (IncrementalMagFreqDist mfd : incrFuncs)
//			if (mfd != crustalFault)
			sum.addIncrementalMagFreqDist(mfd);
		sum.setName("Total");
		incrFuncs.add(0, sum);
		cmlFuncs.add(0, sum.getCumRateDistWithOffset());
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
		texFuncs.add(cmlFuncs.get(0));
		texPrefixes.add("Total");
		
		SummedMagFreqDist gridSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (IncrementalMagFreqDist mfd : incrFuncs)
			if (mfd.getName().contains("Gridded") || mfd.getName().contains("Slab"))
			gridSum.addIncrementalMagFreqDist(mfd);
		gridSum.setName("Total Gridded");
		incrFuncs.add(1, gridSum);
		cmlFuncs.add(1, gridSum.getCumRateDistWithOffset());
		chars.add(1, new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLACK));
		texFuncs.add(cmlFuncs.get(1));
		texPrefixes.add("TotalGridded");
		
		for (boolean includeObs : new boolean[] {false,true}) {
			List<IncrementalMagFreqDist> myIncrFuncs = new ArrayList<>();
			List<DiscretizedFunc> myCmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> myChars = new ArrayList<>();
			
			if (includeObs) {
				boolean useRateModelUncert = true;
				boolean use1973 = false;
				
				double obsMmax = 0d;
				int transAlpha = 60;
				
				boolean[] rateBools = use1973 ? new boolean[] {true,false} : new boolean[] {true};
				for (boolean rates1900 : rateBools) {
					InputStream is = SeismicityRateFileLoader.class.getResourceAsStream(
							"/data/erf/prvi25/seismicity/rates/directrates_2025_05_08/directrates-PRVI Union-Full-"
									+(rates1900 ? "1900" : "1973")+"-2024.csv");
					CSVFile<String> totalRateCSV = CSVFile.readStream(is, false);
					List<Direct> directs = SeismicityRateFileLoader.loadDirectBranches(totalRateCSV);

					IncrementalMagFreqDist obsRefMFD;
					if (rates1900)
//						obsRefMFD = FaultSysTools.initEmptyMFD(6.01, directs.get(0).incrementalDist.getMaxX()-0.1);
						obsRefMFD = FaultSysTools.initEmptyMFD(6.01, directs.get(0).maxObsIncrMag);
					else
						obsRefMFD = FaultSysTools.initEmptyMFD(5.01, 5.99);
					EvenlyDiscretizedFunc refCml = obsRefMFD.getCumRateDistWithOffset();
					
					obsMmax = Math.max(obsMmax, obsRefMFD.getMaxX());
					
					IncrementalMagFreqDist meanIncrObs = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
					IncrementalMagFreqDist incr2p5 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
					IncrementalMagFreqDist incr16 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
					IncrementalMagFreqDist incr84 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
					IncrementalMagFreqDist incr97p5 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
					
					IncrementalMagFreqDist[] dataIncrs = {meanIncrObs, incr2p5, incr16, incr84, incr97p5};
					
					EvenlyDiscretizedFunc meanObsCml = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
					EvenlyDiscretizedFunc cml2p5 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
					EvenlyDiscretizedFunc cml16 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
					EvenlyDiscretizedFunc cml84 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
					EvenlyDiscretizedFunc cml97p5 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
					
					EvenlyDiscretizedFunc[] dataCmls = {meanObsCml, cml2p5, cml16, cml84, cml97p5};
					
					for (int i=0; i<meanIncrObs.size(); i++) {
						double mag = meanIncrObs.getX(i);
						Preconditions.checkState((float)mag >= (float)directs.get(0).M1, "mag=%s, M1=%s", (float)mag, (float)directs.get(0).M1);
						if ((float)mag > (float)directs.get(0).incrementalDist.getMaxX()) {
							break;
						}
						Preconditions.checkState(directs.size() == dataIncrs.length);
						for (int j=0; j<directs.size(); j++) {
							EvenlyDiscretizedFunc directIncr = directs.get(j).incrementalDist;
							EvenlyDiscretizedFunc directCml = directs.get(j).cumulativeDist;
							
							dataIncrs[j].set(i, directIncr.getY(directIncr.getClosestXIndex(mag)));
							dataCmls[j].set(i, directCml.getY(directCml.getClosestXIndex(mag-0.05)));
						}
					}
					
					meanIncrObs.setName(rates1900 ? "Observed" : null);
					myIncrFuncs.add(meanIncrObs);
					meanObsCml.setName(meanIncrObs.getName());
					myCmlFuncs.add(meanObsCml);
					myChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY));
					
					if (!useRateModelUncert) {
						UncertainBoundedIncrMagFreqDist incrBounds95 = new UncertainBoundedIncrMagFreqDist(
								meanIncrObs, incr2p5, incr97p5, UncertaintyBoundType.CONF_95);
						UncertainBoundedIncrMagFreqDist incrBounds68 = new UncertainBoundedIncrMagFreqDist(
								meanIncrObs, incr16, incr84, UncertaintyBoundType.CONF_68);
						UncertainArbDiscFunc cmlBounds95 = new UncertainArbDiscFunc(meanObsCml, cml2p5, cml97p5);
						UncertainArbDiscFunc cmlBounds68 = new UncertainArbDiscFunc(meanObsCml, cml16, cml84);
						
						incrBounds95.setName(rates1900 ? "68% and 95% bounds" : null);
						myIncrFuncs.add(incrBounds95);
						
						cmlBounds95.setName(incrBounds95.getName());
						myCmlFuncs.add(cmlBounds95);
						myChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 2f, new Color(0, 0, 0, transAlpha)));
						
						incrBounds68.setName(null);
						myIncrFuncs.add(incrBounds68);
						cmlBounds68.setName(null);
						myCmlFuncs.add(cmlBounds68);
						myChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 2f, new Color(0, 0, 0, transAlpha)));
					}
				}
				
				if (useRateModelUncert) {
					InputStream is = SeismicityRateFileLoader.class.getResourceAsStream(
//							"/data/erf/prvi25/seismicity/rates/directrates_2025_05_08/rateunc-Union-Full-v3.csv");
							"/data/erf/prvi25/seismicity/rates/directrates_2025_05_08/rateunc-mix-Union-Full-v3.csv");
					CSVFile<String> unionRates = CSVFile.readStream(is, false);
					List<Exact> uncertBranches = SeismicityRateFileLoader.loadExactBranches(unionRates);
					
					EvenlyDiscretizedFunc cml2p5 = SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.025).cumulativeDist;
					EvenlyDiscretizedFunc cml16 = SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.16).cumulativeDist;
					EvenlyDiscretizedFunc cml84 = SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.84).cumulativeDist;
					EvenlyDiscretizedFunc cml97p5 = SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.975).cumulativeDist;
					
					IncrementalMagFreqDist obsRefMFD = FaultSysTools.initEmptyMFD(5.01, obsMmax);
					IncrementalMagFreqDist incr2p5 = SeismicityRateFileLoader.buildIncrementalMFD(
							SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.025), obsRefMFD, obsMmax); 
					IncrementalMagFreqDist incr16 = SeismicityRateFileLoader.buildIncrementalMFD(
							SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.16), obsRefMFD, obsMmax);
					IncrementalMagFreqDist incr84 = SeismicityRateFileLoader.buildIncrementalMFD(
							SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.84), obsRefMFD, obsMmax);
					IncrementalMagFreqDist incr97p5 = SeismicityRateFileLoader.buildIncrementalMFD(
							SeismicityRateFileLoader.locateQuantile(uncertBranches, 0.975), obsRefMFD, obsMmax);
					
					IncrementalMagFreqDist averageIncr = average(incr16, incr84);
					EvenlyDiscretizedFunc averageCml = averageCml(cml16, cml84);
					
					UncertainBoundedIncrMagFreqDist incrBounds95 = new UncertainBoundedIncrMagFreqDist(
							averageIncr, incr2p5, incr97p5, UncertaintyBoundType.CONF_95);
					UncertainBoundedIncrMagFreqDist incrBounds68 = new UncertainBoundedIncrMagFreqDist(
							averageIncr, incr16, incr84, UncertaintyBoundType.CONF_68);
					UncertainArbDiscFunc cmlBounds95 = new UncertainArbDiscFunc(averageCml, cml2p5, cml97p5);
					UncertainArbDiscFunc cmlBounds68 = new UncertainArbDiscFunc(averageCml, cml16, cml84);
					
					incrBounds95.setName("68% and 95% bounds");
					myIncrFuncs.add(incrBounds95);
					cmlBounds95.setName(incrBounds95.getName());
					myCmlFuncs.add(cmlBounds95);
					myChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 2f, new Color(0, 0, 0, transAlpha)));
					
					incrBounds68.setName(null);
					myIncrFuncs.add(incrBounds68);
					cmlBounds68.setName(null);
					myCmlFuncs.add(cmlBounds68);
					myChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 2f, new Color(0, 0, 0, transAlpha)));
				}
			}
			
			myIncrFuncs.addAll(incrFuncs);
			myCmlFuncs.addAll(cmlFuncs);
			myChars.addAll(chars);
			
			PlotSpec incrSpec = new PlotSpec(myIncrFuncs, myChars, " ", "Magnitude", "Incremental Rate (1/yr)");
			incrSpec.setLegendInset(includeObs ? RectangleAnchor.BOTTOM_LEFT : RectangleAnchor.TOP_RIGHT);
			PlotSpec cmlSpec = new PlotSpec(myCmlFuncs, myChars, " ", "Magnitude", "Cumulative Rate (1/yr)");
			cmlSpec.setLegendInset(includeObs ? RectangleAnchor.BOTTOM_LEFT : RectangleAnchor.TOP_RIGHT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			if (includeObs)
				gp.getPlotPrefs().setLegendFontSize(16);

			Range yRange = new Range(1e-6, 1e1);
			Range xRange = new Range(5d, 9.5d);
			
			String prefix = "combined_mfds";
			if (includeObs)
				prefix += "_obs";

			gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
			gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
			PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, true, false);
		}
		
		Preconditions.checkState(texFuncs.size() == texPrefixes.size());
		
		double[] texMags = {5d, 6d, 7d, 8d, 9d};
		String[] texMagLabels = { "Five", "Six", "Seven", "Eight", "Nine" };
		FileWriter texFW = new FileWriter(new File(outputDir, "combined_mfds.tex"));
		for (int m=0; m<texMags.length; m++) {
			String ratePrefix = "CmlRateM"+texMagLabels[m];
			String riPrefix = "CmlRIM"+texMagLabels[m];
			for (int i=0; i<texFuncs.size(); i++) {
				EvenlyDiscretizedFunc func = texFuncs.get(i);
				String label = texPrefixes.get(i);
				double rate = func.getY(texMags[m]);
				double ri = 1d/rate;
				System.out.println(label+" M>"+(float)texMags[m]+": "+(float)rate+" ("+(float)ri+")");
				if (rate == 0d) {
//					texFW.write(LaTeXUtils.defineValueCommand(ratePrefix+label, "0")+"\n");
//					texFW.write(LaTeXUtils.defineValueCommand(riPrefix+label, "$\\infty$", false)+"\n");
//					texFW.write(LaTeXUtils.defineValueCommand(ratePrefix+label, "")+"\n");
//					texFW.write(LaTeXUtils.defineValueCommand(riPrefix+label, "")+"\n");
				} else {
					texFW.write(LaTeXUtils.defineValueCommand(ratePrefix+label,
							LaTeXUtils.numberExpFormatSigFigs(rate, 2), false)+"\n");
					if (ri > 5)
						texFW.write(LaTeXUtils.defineValueCommand(riPrefix+label,
								LaTeXUtils.groupedIntNumber(ri), false)+"\n");
					else
						texFW.write(LaTeXUtils.defineValueCommand(riPrefix+label,
								LaTeXUtils.numberExpFormatFixedDecimal(ri, 1), false)+"\n");
				}
			}
		}
		texFW.close();
	}
	
	static IncrementalMagFreqDist calcFaultMFD(Region region, FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) {
		return sol.calcNucleationMFD_forRegion(region, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
	}
	
	static IncrementalMagFreqDist calcGriddedMFD(Region region, TectonicRegionType trt,
			FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) {
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		SummedMagFreqDist ret = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int i=0; i<gridProv.getNumLocations(); i++) {
			if (region.contains(gridProv.getLocation(i))) {
				IncrementalMagFreqDist mfd = gridProv.getMFD(trt, i);
				ret.addIncrementalMagFreqDist(mfd);
			}
		}
		return ret;
	}
	
	static IncrementalMagFreqDist average(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		Preconditions.checkState(mfd2.size() == mfd1.size());
		Preconditions.checkState(mfd2.getMinX() == mfd1.getMinX());
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(mfd1.getMinX(), mfd1.size(), mfd1.getDelta());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, 0.5*mfd1.getY(i) + 0.5*mfd2.getY(i));
		return ret;
	}
	
	static EvenlyDiscretizedFunc averageCml(EvenlyDiscretizedFunc mfd1, EvenlyDiscretizedFunc mfd2) {
		Preconditions.checkState(mfd2.size() == mfd1.size());
		Preconditions.checkState(mfd2.getMinX() == mfd1.getMinX());
		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(mfd1.getMinX(), mfd1.size(), mfd1.getDelta());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, 0.5*mfd1.getY(i) + 0.5*mfd2.getY(i));
		return ret;
	}
	
	static IncrementalMagFreqDist sum(IncrementalMagFreqDist... mfds) {
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(mfds[0].getMinX(), mfds[0].size(), mfds[0].getDelta());
		for (IncrementalMagFreqDist mfd : mfds) {
			Preconditions.checkState(ret.size() == mfd.size());
			Preconditions.checkState(ret.getMinX() == mfd.getMinX());
			for (int i=0; i<ret.size(); i++)
				ret.add(i, mfd.getY(i));
		}
		return ret;
	}

}
