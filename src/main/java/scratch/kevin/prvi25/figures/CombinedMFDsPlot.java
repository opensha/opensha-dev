package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Direct;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Exact;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
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
//		FaultSystemSolution subductionSol1 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2025_05_21-prvi25_crustal_subduction_combined_branches-ba_only-slab_mc_7p4-vs760/combined_branch_averaged_solution.zip"));
		FaultSystemSolution subductionSol2 = FaultSystemSolution.load(SUBDUCTION_SOL_SMALL);
		
		File outputDir = new File(FIGURES_DIR, "combined_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		boolean useRateModelUncert = true;
//		boolean use1973 = false;
		boolean plot1973 = true;
		boolean plot1900 = true;
		boolean plotEpochOverlap = true;
		boolean includeTotalGridded = false;
		
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(9.45);
		
		List<PRVI25_SeismicityRegions> seisRegions = new ArrayList<>(List.of(PRVI25_SeismicityRegions.values()));
		seisRegions.add(null); // total
		
		List<EvenlyDiscretizedFunc> texFuncs = new ArrayList<>();
		List<String> texPrefixes = new ArrayList<>();
		
		IncrementalMagFreqDist crustalFault = calcFaultMFD(PRVI25_SeismicityRegions.CRUSTAL.load(), crustalSol, refMFD);
		EvenlyDiscretizedFunc crustalFaultCml = crustalFault.getCumRateDistWithOffset();
		texFuncs.add(crustalFaultCml);
		texPrefixes.add("CrustalFault");
		
		IncrementalMagFreqDist crustalGrid = calcGriddedMFD(PRVI25_SeismicityRegions.CRUSTAL, TectonicRegionType.ACTIVE_SHALLOW, crustalSol, refMFD);
		EvenlyDiscretizedFunc crustalGridCml = crustalGrid.getCumRateDistWithOffset();
		texFuncs.add(crustalGridCml);
		texPrefixes.add("CrustalGridded");
		
		IncrementalMagFreqDist crustalTotal = sum(crustalFault, crustalGrid);
		EvenlyDiscretizedFunc crustalTotalCml = crustalTotal.getCumRateDistWithOffset();
		texFuncs.add(crustalTotalCml);
		texPrefixes.add("CrustalTotal");
		
//		Region mueInterfaceReg = PRVI25_SeismicityRegions.MUE_INTERFACE.load();
		// bigger polygon to make sure we don't accidentally reduce any rates for section points barely outside of seis region
		Region mueInterfaceReg = new Region(LocationList.of(
				new Location(19, -71),
				new Location(18, -64),
				new Location(16, -64),
				new Location(16, -71),
				new Location(19, -71)),
				BorderType.MERCATOR_LINEAR);
		IncrementalMagFreqDist faultMuertos = average(
				calcFaultMFD(mueInterfaceReg, subductionSol1, refMFD),
				calcFaultMFD(mueInterfaceReg, subductionSol2, refMFD));
		EvenlyDiscretizedFunc faultMuertosCml = faultMuertos.getCumRateDistWithOffset();
		texFuncs.add(faultMuertosCml);
		texPrefixes.add("MuertosInterfaceFault");
		
		IncrementalMagFreqDist griddedMuertos = average(
				calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTERFACE, TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol1, refMFD),
				calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTERFACE, TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol2, refMFD));
		EvenlyDiscretizedFunc griddedMuertosCml = griddedMuertos.getCumRateDistWithOffset();
		texFuncs.add(griddedMuertosCml);
		texPrefixes.add("MuertosInterfaceGridded");
		
		IncrementalMagFreqDist mueTotal = sum(faultMuertos, griddedMuertos);
		EvenlyDiscretizedFunc mueTotalCml = mueTotal.getCumRateDistWithOffset();
		texFuncs.add(mueTotalCml);
		texPrefixes.add("MuertosInterfaceTotal");
		
		IncrementalMagFreqDist slabMuertos = calcGriddedMFD(PRVI25_SeismicityRegions.MUE_INTRASLAB,
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol1, refMFD);
		EvenlyDiscretizedFunc slabMuertosCml = slabMuertos.getCumRateDistWithOffset();
		texFuncs.add(slabMuertosCml);
		texPrefixes.add("MuertosSlab");
		
//		Region carInterfaceReg = CAR_INTERFACE.MUE_INTERFACE.load();
		// bigger polygon to make sure we don't accidentally reduce any rates for section points barely outside of seis region
		Region carInterfaceReg = new Region(LocationList.of(
				new Location(21, -71),
				new Location(21, -61.7),
				new Location(16, -61.7),
				new Location(18, -64),
				new Location(19, -71),
				new Location(21, -71)),
				BorderType.MERCATOR_LINEAR);
		IncrementalMagFreqDist faultCar = average(
				calcFaultMFD(carInterfaceReg, subductionSol1, refMFD),
				calcFaultMFD(carInterfaceReg, subductionSol2, refMFD));
		EvenlyDiscretizedFunc faultCarCml = faultCar.getCumRateDistWithOffset();
		texFuncs.add(faultCarCml);
		texPrefixes.add("CarInterfaceFault");
		
		IncrementalMagFreqDist griddedCar = average(
				calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTERFACE, TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol1, refMFD),
				calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTERFACE, TectonicRegionType.SUBDUCTION_INTERFACE, subductionSol2, refMFD));
		EvenlyDiscretizedFunc griddedCarCml = griddedCar.getCumRateDistWithOffset();
		texFuncs.add(griddedCarCml);
		texPrefixes.add("CarInterfaceGridded");
		
		IncrementalMagFreqDist carTotal = sum(faultCar, griddedCar);
		EvenlyDiscretizedFunc carTotalCml = carTotal.getCumRateDistWithOffset();
		texFuncs.add(carTotalCml);
		texPrefixes.add("CarInterfaceTotal");
		
		IncrementalMagFreqDist slabCar = calcGriddedMFD(PRVI25_SeismicityRegions.CAR_INTRASLAB,
				TectonicRegionType.SUBDUCTION_SLAB, subductionSol1, refMFD);
		EvenlyDiscretizedFunc slabCarCml = slabCar.getCumRateDistWithOffset();
		texFuncs.add(slabCarCml);
		texPrefixes.add("CarSlab");
		
		texFuncs.add(sum(faultCar, griddedCar, faultMuertos, griddedMuertos).getCumRateDistWithOffset());
		texPrefixes.add("InterfaceTotal");
		
		IncrementalMagFreqDist slabTotal = sum(slabCar, slabMuertos);
		texFuncs.add(slabTotal.getCumRateDistWithOffset());
		texPrefixes.add("SlabTotal");
		
		GridSourceList gridProv = subductionSol1.requireModule(GridSourceList.class);
		SummedMagFreqDist slabSumTest = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int i=0; i<gridProv.getNumLocations(); i++) {
			IncrementalMagFreqDist nodeMFD = gridProv.getMFD(TectonicRegionType.SUBDUCTION_SLAB, i);
			if (nodeMFD != null)
				slabSumTest.addIncrementalMagFreqDist(nodeMFD);
		}
		System.out.println("Sum slab total:\t"+(float)slabTotal.calcSumOfY_Vals());
		System.out.println("Sum slab raw test:\t"+(float)slabSumTest.calcSumOfY_Vals());
		for (int i=0; i<refMFD.size(); i++) {
			double mag = refMFD.getX(i);
			double sum = slabTotal.getY(i);
			double testSum = slabSumTest.getY(i);
			Preconditions.checkState(Precision.equals(sum, testSum, 1e-4),
					"Slab sum test mismatch for M=%s. sum:\t%s;\trawTestSum=%s;\tdiff=%s",
					(float)mag, (float)sum, (float)testSum, (float)(sum-testSum));
		}
		
		for (PRVI25_SeismicityRegions seisReg : seisRegions) {
			String prefix = "combined_mfds";
			if (seisReg != null)
				prefix += "_"+seisReg.name();
			
			List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
			List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			EvenlyDiscretizedFunc cmlSumFunc = null;
			
			if (seisReg == null || seisReg == PRVI25_SeismicityRegions.CRUSTAL) {
				incrFuncs.add(crustalFault);
				cmlFuncs.add(crustalFaultCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue));
				
				incrFuncs.add(crustalGrid);
				cmlFuncs.add(crustalGridCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightblue));
				
				if (seisReg != null) {
					crustalFault.setName("Model (on-fault)");
					crustalGrid.setName("Model (gridded)");
					crustalTotal.setName("Total model");
					
					incrFuncs.add(crustalTotal);
					cmlFuncs.add(crustalTotalCml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					cmlSumFunc = crustalTotalCml;
				} else {
					crustalFault.setName("Crustal, on-fault");
					crustalGrid.setName("Crustal, gridded");
				}
			}
			
			if (seisReg == null || seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE) {
				incrFuncs.add(faultMuertos);
				cmlFuncs.add(faultMuertosCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.green));

				incrFuncs.add(griddedMuertos);
				cmlFuncs.add(griddedMuertosCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightgreen));
				
				if (seisReg != null) {
					faultMuertos.setName("Model (trench-breaking)");
					griddedMuertos.setName("Model (gridded)");
					mueTotal.setName("Total model");
					
					incrFuncs.add(mueTotal);
					cmlFuncs.add(mueTotalCml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					cmlSumFunc = mueTotalCml;
				} else {
					faultMuertos.setName("Muertos, interface trench-breaking");
					griddedMuertos.setName("Muertos, interface gridded");
				}
			}
			
			if (seisReg == null || seisReg == PRVI25_SeismicityRegions.MUE_INTRASLAB) {
				incrFuncs.add(slabMuertos);
				cmlFuncs.add(slabMuertosCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightgreen));
				cmlSumFunc = slabMuertosCml;
				
				if (seisReg != null)
					slabMuertos.setName("Model");
				else
					slabMuertos.setName("Muertos, intraslab");
			}
			
			if (seisReg == null || seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE) {
				incrFuncs.add(faultCar);
				cmlFuncs.add(faultCarCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange));

				incrFuncs.add(griddedCar);
				cmlFuncs.add(griddedCarCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_lightorange));
				
				if (seisReg != null) {
					faultCar.setName("Model (trench-breaking)");
					griddedCar.setName("Model (gridded)");
					carTotal.setName("Total model");
					
					incrFuncs.add(carTotal);
					cmlFuncs.add(carTotalCml);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					cmlSumFunc = carTotalCml;
				} else {
					faultCar.setName("Caribbean, interface trench-breaking");
					griddedCar.setName("Caribbean, interface gridded");
				}
			}
			
			if (seisReg == null || seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB) {
				incrFuncs.add(slabCar);
				cmlFuncs.add(slabCarCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Colors.tab_lightorange));
				cmlSumFunc = slabCarCml;
				
				if (seisReg != null)
					slabCar.setName("Model");
				else
					slabCar.setName("Caribbean, intraslab");
			}
			
			// set cumulatives to same name as incrementals
			for (int i=0; i<incrFuncs.size(); i++)
				cmlFuncs.get(i).setName(incrFuncs.get(i).getName());
			
			boolean[] obsBools;
			
			if (seisReg == null) {
				SummedMagFreqDist sum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				for (IncrementalMagFreqDist mfd : incrFuncs)
//					if (mfd != crustalFault)
					sum.addIncrementalMagFreqDist(mfd);
				sum.setName("Total");
				cmlSumFunc = sum.getCumRateDistWithOffset();
				incrFuncs.add(0, sum);
				cmlFuncs.add(0, cmlSumFunc);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
				texFuncs.add(cmlFuncs.get(0));
				texPrefixes.add("Total");
				
				SummedMagFreqDist gridSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				for (IncrementalMagFreqDist mfd : incrFuncs)
					if (mfd.getName().contains("Gridded") || mfd.getName().contains("Slab"))
					gridSum.addIncrementalMagFreqDist(mfd);
				gridSum.setName("Total Gridded");
				if (includeTotalGridded) {
					incrFuncs.add(1, gridSum);
					cmlFuncs.add(1, gridSum.getCumRateDistWithOffset());
					chars.add(1, new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLACK));
				}
				texFuncs.add(cmlFuncs.get(1));
				texPrefixes.add("TotalGridded");
				
				obsBools = new boolean[] {false,true};
			} else {
				obsBools = new boolean[] {true};
			}
			
			String regName;
			String title;
			String obsTexPrefix;
			if (seisReg == null) {
				regName = "Union";
				title = "Combined model";
				obsTexPrefix = "TotalObs";
			} else if (seisReg == PRVI25_SeismicityRegions.CRUSTAL) {
				regName = "Crustal";
				title = "Crustal";
				obsTexPrefix = "CrustalObs";
			} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE) {
				regName = "CAR Interface";
				title = "Caribbean interface";
				obsTexPrefix = "CarInterfaceObs";
			} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB) {
				regName = "CAR Intraslab";
				title = "Caribbean intraslab";
				obsTexPrefix = "CarSlabObs";
			} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE) {
				regName = "MUE Interface";
				title = "Muertos interface";
				obsTexPrefix = "MuertosInterfaceObs";
			} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTRASLAB) {
				regName = "MUE Intraslab";
				title = "Muertos intraslab";
				obsTexPrefix = "MuertosSlabObs";
			} else {
				throw new IllegalStateException("Unknown seis region: " + seisReg);
			}
			
			for (boolean includeObs : obsBools) {
				List<IncrementalMagFreqDist> myIncrFuncs = new ArrayList<>();
				List<DiscretizedFunc> myCmlFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> myChars = new ArrayList<>();
				
				double obsN5 = Double.NaN;
				double obsN6 = Double.NaN;
				
				if (includeObs) {
					double obsMmax = 0d;
					int transAlpha = 60;
					
					String ratesPrefix = "/data/erf/prvi25/seismicity/rates/directrates_2025_05_08/";
					
					Preconditions.checkState(plot1900 || plot1973);
					boolean[] rateBools;
					if (plot1900 && plot1973)
						rateBools = new boolean[] {true,false};
					else if (plot1900)
						rateBools = new boolean[] {true};
					else
						rateBools = new boolean[] {false};
					plotEpochOverlap &= rateBools.length>1;
					
					EvenlyDiscretizedFunc combObsCml = refMFD.getCumRateDistWithOffset().deepClone();
					Preconditions.checkState(combObsCml.calcSumOfY_Vals() == 0d);
					for (boolean is1900 : rateBools) {
						String ratesPath = ratesPrefix+"directrates-PRVI "+regName+"-Full-"+(is1900 ? "1900" : "1973")+"-2024.csv";
						System.out.println("Loading direct rates from "+ratesPath);
						InputStream is = SeismicityRateFileLoader.class.getResourceAsStream(ratesPath);
						CSVFile<String> totalRateCSV = CSVFile.readStream(is, false);
						List<Direct> directs = SeismicityRateFileLoader.loadDirectBranches(totalRateCSV);
						
						double minFuncMag = is1900 ? 6.01 : 5.01;
						double maxFuncMag = directs.get(0).maxObsIncrMag;
						if (maxFuncMag == 0d) {
							// not found
//							maxFuncMag = directs.get(0).cumulativeDist.getMaxX()+0.01;
							maxFuncMag = directs.get(0).cumulativeDist.getMaxX()-0.01;
						}
						if (!is1900 && plot1900 && !plotEpochOverlap)
							// if these are 1973 rates and we're including 1900 and we're not plotting overlap, trim before 1900 rates
							maxFuncMag = Math.min(maxFuncMag, 5.99);

						IncrementalMagFreqDist obsRefMFD = FaultSysTools.initEmptyMFD(minFuncMag, maxFuncMag);
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
						
						if (is1900)
							obsN6 = meanObsCml.getY(0);
						else
							obsN5 = meanObsCml.getY(0);
						
						for (int i=0; i<meanObsCml.size(); i++) {
							double x = meanObsCml.getX(i);
							if (plotEpochOverlap && !is1900 && x >= 5.99)
								// we're including 1900 rates, stop these 1973 rates below M6
								break;
							combObsCml.set(combObsCml.getClosestXIndex(x), meanObsCml.getY(i));
						}
						
						PlotCurveCharacterstics obsChar;
						String obsName;
						Color darkishGray = new Color(80, 80, 80);
						if (plotEpochOverlap) {
							// we're including both in the plot, and they'll overlap
							if (is1900) {
								obsName = "Observed (1900-2023)";
								obsChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, Color.BLACK);
							} else {
								obsName = "Observed (1973-2023)";
								obsChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, darkishGray);
							}
						} else if (plot1900 && plot1973) {
							// we're including both as one stitched line
							if (is1900)
								obsName = "Observed";
							else
								obsName = null; // hide label
							obsChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, darkishGray);
						} else {
							if (is1900)
								obsName = "Observed (1900-2023)";
							else
								obsName = "Observed (1973-2023)";
							obsChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, darkishGray);
						}
						meanIncrObs.setName(obsName);
						myIncrFuncs.add(meanIncrObs);
						meanObsCml.setName(meanIncrObs.getName());
						myCmlFuncs.add(meanObsCml);
						myChars.add(obsChar);
						
						if (!useRateModelUncert) {
							UncertainBoundedIncrMagFreqDist incrBounds95 = new UncertainBoundedIncrMagFreqDist(
									meanIncrObs, incr2p5, incr97p5, UncertaintyBoundType.CONF_95);
							UncertainBoundedIncrMagFreqDist incrBounds68 = new UncertainBoundedIncrMagFreqDist(
									meanIncrObs, incr16, incr84, UncertaintyBoundType.CONF_68);
							UncertainArbDiscFunc cmlBounds95 = new UncertainArbDiscFunc(meanObsCml, cml2p5, cml97p5);
							UncertainArbDiscFunc cmlBounds68 = new UncertainArbDiscFunc(meanObsCml, cml16, cml84);
							
							incrBounds95.setName(is1900 ? "68% and 95% bounds" : null);
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
					
					texFuncs.add(combObsCml);
					texPrefixes.add(obsTexPrefix);
					System.out.println(combObsCml);
					
					if (useRateModelUncert) {
						List<Exact> uncertBranches;
						if (seisReg == null) {
							InputStream is = SeismicityRateFileLoader.class.getResourceAsStream(
									ratesPrefix+"rateunc-Union-Full-v3.csv");
//									ratesPrefix+"rateunc-1900-M6-Union-Full-v4.csv");
//									ratesPrefix+"rateunc-mix-Union-Full-v3.csv");
							CSVFile<String> unionRates = CSVFile.readStream(is, false);
							uncertBranches = SeismicityRateFileLoader.loadExactBranches(unionRates);
						} else if (seisReg == PRVI25_SeismicityRegions.CRUSTAL) {
							uncertBranches = (List<Exact>)PRVI25_CrustalSeismicityRate.loadRates(RateType.EXACT);
						} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE) {
							uncertBranches = (List<Exact>)PRVI25_SubductionCaribbeanSeismicityRate.loadRates(RateType.EXACT, false);
						} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB) {
							uncertBranches = (List<Exact>)PRVI25_SubductionCaribbeanSeismicityRate.loadRates(RateType.EXACT, true);
						} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE) {
							uncertBranches = (List<Exact>)PRVI25_SubductionMuertosSeismicityRate.loadRates(RateType.EXACT, false);
						} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTRASLAB) {
							uncertBranches = (List<Exact>)PRVI25_SubductionMuertosSeismicityRate.loadRates(RateType.EXACT, true);
						} else {
							throw new IllegalStateException();
						}
						
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
				
				int obsNum = myIncrFuncs.size();
				
				myIncrFuncs.addAll(incrFuncs);
				myCmlFuncs.addAll(cmlFuncs);
				myChars.addAll(chars);
				
				if (includeObs) {
					for (int i=0; i<obsNum; i++) {
						PlotCurveCharacterstics obsChar = myChars.get(i);
						if (obsChar.getLineType() != PlotLineType.SHADED_UNCERTAIN) {
							// add the obs again on top
							IncrementalMagFreqDist obsCopy = myIncrFuncs.get(i).deepClone();
							obsCopy.setName(null);
							DiscretizedFunc obsCmlCopy = myCmlFuncs.get(i).deepClone();
							obsCmlCopy.setName(null);
							myIncrFuncs.add(obsCopy);
							myCmlFuncs.add(obsCmlCopy);
							myChars.add(obsChar);
						}
					}
				}
				
				PlotSpec incrSpec = new PlotSpec(myIncrFuncs, myChars, " ", "Magnitude", "Incremental Rate (1/yr)");
				incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				PlotSpec cmlSpec = new PlotSpec(myCmlFuncs, myChars, " ", "Magnitude", "Cumulative Rate (1/yr)");
				cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				if (includeObs) {
					Preconditions.checkNotNull(cmlSumFunc);
					double modelN5 = cmlSumFunc.getY(cmlSumFunc.getClosestXIndex(5.01));
					double modelN6 = cmlSumFunc.getY(cmlSumFunc.getClosestXIndex(6.01));
					Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
					Font titleFont = new Font(Font.SANS_SERIF, Font.BOLD, 26);
					
					DecimalFormat df5 = new DecimalFormat("0.00");
					DecimalFormat df6 = new DecimalFormat("0.00");
					
					String modelText = "Model:  M≥5="+df5.format(modelN5)+";  M≥6="+df6.format(modelN6);
					String obsText = "Observed:  M≥5="+df5.format(obsN5)+";  M≥6="+df6.format(obsN6);

					XYTextAnnotation regAnn = new XYTextAnnotation(title, 9.3, 7e0);
					regAnn.setFont(titleFont);
					regAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
					cmlSpec.addPlotAnnotation(regAnn);
					XYTextAnnotation modelAnn = new XYTextAnnotation(modelText, 9.3, 2.45e0);
					modelAnn.setFont(annFont);
					modelAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
					cmlSpec.addPlotAnnotation(modelAnn);
					XYTextAnnotation obsAnn = new XYTextAnnotation(obsText, 9.3, 1.1e0);
					obsAnn.setFont(annFont);
					obsAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
					cmlSpec.addPlotAnnotation(obsAnn);
				}
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				if (includeObs)
					gp.getPlotPrefs().setLegendFontSize(16);
				
//				if (seisReg != null)
//					gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);

				Range yRange = new Range(1e-6, 1e1);
				Range xRange = new Range(5d, 9.5d);
				
				String myPrefix = prefix;
				if (includeObs)
					myPrefix += "_obs";
				
				int smallW = 600;
				int smallH = 550;
				int largeW = 800;
				int largeH = 750;
				
				if (seisReg != null) {
					// only do small cumulative plots for region-specific
					gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
					PlotUtils.writePlots(outputDir, myPrefix+"_cml", gp, smallW, smallH, true, true, false);
				} else {
					// do incremental (large)
					gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
					PlotUtils.writePlots(outputDir, myPrefix, gp, largeW, largeH, true, true, false);
					// both small and large for total
					gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
					PlotUtils.writePlots(outputDir, myPrefix+"_cml", gp, largeW, largeH, true, true, false);
					// clear out all labels except for total and data
					// start at obsNum+1 to keep total
					for (int i=obsNum+1; i<myCmlFuncs.size(); i++) {
						DiscretizedFunc func = myCmlFuncs.get(i).deepClone();
						func.setName(null);
						myCmlFuncs.set(i, func);
					}
					gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
					PlotUtils.writePlots(outputDir, myPrefix+"_cml_small", gp, smallW, smallH, true, true, false);
				}
			}
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
				Preconditions.checkState(!(func instanceof IncrementalMagFreqDist),
						"TexFunc %s (%s) is an incr MFD", i, func.getName());
				String label = texPrefixes.get(i);
				int magIndex = func.getClosestXIndex(texMags[m]);
				Preconditions.checkState((float)func.getX(magIndex) == (float)texMags[m], "Bad mag compare %s != %s", (float)func.getX(magIndex), (float)texMags[m]);
				double rate = func.getY(magIndex);
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
	
	static IncrementalMagFreqDist calcGriddedMFD(PRVI25_SeismicityRegions seisReg, TectonicRegionType trt,
			FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) throws IOException {
		GriddedGeoDataSet depths = null;
		GriddedGeoDataSet altDepths = null;
		GriddedRegion gridReg;
		if (trt == TectonicRegionType.SUBDUCTION_SLAB) {
			depths = PRVI25_GridSourceBuilder.loadSubductionDepths(seisReg);
			// store the other one in alt
			if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB)
				altDepths = PRVI25_GridSourceBuilder.loadSubductionDepths(PRVI25_SeismicityRegions.MUE_INTRASLAB);
			else
				altDepths = PRVI25_GridSourceBuilder.loadSubductionDepths(PRVI25_SeismicityRegions.CAR_INTRASLAB);
			gridReg = depths.getRegion();
		} else {
			Region region = seisReg.load();
			gridReg = new GriddedRegion(region, 0.1, GriddedRegion.ANCHOR_0_0);
		}
		double maxSlabDepthMismatch = 0d;
		GridSourceList gridProv = sol.requireModule(GridSourceList.class);
		SummedMagFreqDist ret = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int i=0; i<gridProv.getNumLocations(); i++) {
			Location loc = gridProv.getLocation(i);
			int regIndex = gridReg.indexForLocation(loc);
			if (regIndex >= 0) {
				if (trt == TectonicRegionType.SUBDUCTION_SLAB) {
					// need to compare depths
					HashSet<Float> uniqueDepths = new HashSet<>(2);
					double modelDepth = depths.get(regIndex);
					int altRegIndex = altDepths.indexOf(loc);
					double altModelDepth = altRegIndex >= 0 ? altDepths.get(altRegIndex) : Double.NaN;
//					double depthTol = 1e-4;
					double depthTol = 1d;
					Preconditions.checkState(!Precision.equals(modelDepth, altModelDepth, depthTol),
									"Model and alt model depths collide for %s. %s == %s",
									(float)modelDepth, (float)altModelDepth);
					for (GriddedRupture rup : gridProv.getRuptures(trt, i)) {
						uniqueDepths.add((float)rup.properties.upperDepth);
						Preconditions.checkState(Precision.equals(rup.properties.upperDepth, rup.properties.lowerDepth),
								"Slab upper != lower: %s != %s", (float)rup.properties.upperDepth, (float)rup.properties.lowerDepth);
						if (Precision.equals(rup.properties.upperDepth, modelDepth, depthTol)) {
							ret.add(refMFD.getClosestXIndex(rup.properties.magnitude), rup.rate);
							maxSlabDepthMismatch = Math.max(maxSlabDepthMismatch, Math.abs(rup.properties.upperDepth-modelDepth));
						} else {
							Preconditions.checkState(Precision.equals(rup.properties.upperDepth, altModelDepth, depthTol),
									"Unexpected slab depth at %s. %s is %s, alt is %s, we have %s",
									loc, seisReg.name(), (float)modelDepth, (float)altModelDepth, rup.properties.upperDepth);
						}
					}
					Preconditions.checkState(uniqueDepths.size() <= 2,
							"Expected 1 or 2 depths for %s, have %s, model=%s, values=%s", loc, uniqueDepths.size(), modelDepth, uniqueDepths);
				} else {
					// no overlap, just use MFD method
					IncrementalMagFreqDist mfd = gridProv.getMFD(trt, i);
					ret.addIncrementalMagFreqDist(mfd);
				}
			}
		}
		if (maxSlabDepthMismatch > 0d)
			System.out.println("Max slab depth mismatch was "+(float)maxSlabDepthMismatch);
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
