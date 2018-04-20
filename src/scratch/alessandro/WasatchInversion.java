package scratch.alessandro;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.SegRateConstraint;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

/**
 * This class reads Wasatch inversion data from files, provides methods for getting the various constraints, and runs the inversion
 * @author field n.
 *
 */
public class WasatchInversion {

	final static boolean D = false;	// debugging flag
	
	public final static String ROOT_PATH = "src/scratch/alessandro/";
	final static String ROOT_DATA_DIR = "src/scratch/alessandro/data/"; // where to find the data

	
	// These values are the same for all fault sections
	final static double UPPER_SEIS_DEPTH = 0;
	final static double LOWER_SEIS_DEPTH = 15;
	final static double FAULT_DIP = 45;
	final static double FAULT_RAKE = -90;
	
	
	ArrayList<FaultSectionPrefData> faultSectionDataList;
	ArrayList<SegRateConstraint> sectionRateConstraints;
	int[][] rupSectionMatrix;
	
	final static String SLIP_RATE_FILENAME = "sliprate_wasatch_mod.txt";
	final static String PALEO_RATE_FILENAME = "paleorate_wasatch_and_95p.txt";
	final static String FAULT_TRACE_DIR_NAME = "subsections_traces/";

	public WasatchInversion() {
		
		readData();
		
	}
		
		
	private void readData() {
	
		int numSections, numRuptures;
		double[] sectSlipRate, sectSlipRateStdDev;

		// Read Data:

		try {
			File file = new File(ROOT_DATA_DIR+SLIP_RATE_FILENAME);
			List<String> fileLines = Files.readLines(file, Charset.defaultCharset());
			numSections = fileLines.size();
			if(D) System.out.println("numSections="+numSections);
			int sectIndex = 0;
			sectSlipRate = new double[numSections];
			sectSlipRateStdDev = new double[numSections];
			for (String line : fileLines) {
				//			System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 3, "Expected 3 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				int testIndex = Integer.valueOf(split[0]);
				if(sectIndex != testIndex)
					throw new RuntimeException("Bad section index; "+sectIndex+" != "+testIndex);
				sectSlipRate[sectIndex] = Double.valueOf(split[1]);	
				sectSlipRateStdDev[sectIndex] = Double.valueOf(split[2]);	
				if(D) System.out.println(sectIndex+"\t"+sectSlipRate[sectIndex]+"\t"+sectSlipRateStdDev[sectIndex]);
				sectIndex+=1;
			}


			// Now read section rate constraints
			sectionRateConstraints   = new ArrayList<SegRateConstraint>();
			file = new File(ROOT_DATA_DIR+PALEO_RATE_FILENAME);
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				//			System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length == 5, "Expected 3 items, got %s", split.length);
				//			System.out.println(split[0]+"\t"+split[1]+"\t"+split[2]);
				sectIndex = Integer.valueOf(split[0]);
				double meanRate = Double.valueOf(split[1]);
				double stdDev = Double.valueOf(split[2]);
				double upp95 = Double.valueOf(split[3]);
				double low95 = Double.valueOf(split[4]);
				stdDev = (upp95-low95)/4.0;  // temporary fix until the stdDev values are fixed
				SegRateConstraint sectionRateConstraint = new SegRateConstraint("Section "+split[0]); // Names are not unique!
				sectionRateConstraint.setSegRate(sectIndex, meanRate, stdDev, low95, upp95);
				sectionRateConstraints.add(sectionRateConstraint);
				if (D) System.out.println(sectionRateConstraint.getFaultName()+"\t"+
						sectionRateConstraint.getSegIndex()+"\t"+
						sectionRateConstraint.getMean()+"\t"+
						sectionRateConstraint.getStdDevOfMean()+"\t"+
						sectionRateConstraint.getLower95Conf()+"\t"+
						sectionRateConstraint.getUpper95Conf());
			}


			// Now read section rupture matrix
			file = new File("src/scratch/alessandro/data/Gsr_matrix_new.txt");
			fileLines = Files.readLines(file, Charset.defaultCharset());
			numRuptures = fileLines.size()-1;
			rupSectionMatrix = new int[numSections][numRuptures];
			int rupIndex=-1;
			for (String line : fileLines) {
				if(rupIndex==-1) {
					rupIndex+=1;
					continue;
				}
				//					System.out.println(line);
				line = line.trim();
				String[] split = line.split("\t");	// tab delimited
				Preconditions.checkState(split.length-1 == numSections, "Number of columns (%s) not consistent with numSections (%s)", split.length, numSections);

				int testRupIndex = Integer.valueOf(split[0]);
				Preconditions.checkState(testRupIndex == rupIndex, "Rup index problem on input file (%s vs %s)", testRupIndex, rupIndex);

				for(int s=0; s<numSections; s++) {
					int sectInRup = Integer.valueOf(split[s+1]);
					if(sectInRup == 1) {
						rupSectionMatrix[s][rupIndex] = 1;
					}
				}

				if(D) {
					System.out.println("Rupture "+rupIndex);
					for(int s=0; s<numSections; s++)
						System.out.print(rupSectionMatrix[s][rupIndex]+" ");
					System.out.print("\n");
				}

				rupIndex+=1;
			}


			// Now read traces and make FaultSectionPrefData for each section
			
			faultSectionDataList = new ArrayList<FaultSectionPrefData>();
			for(int s=0; s<numSections; s++) {
				String traceFileName;
				FaultTrace fltTrace = new FaultTrace("Trace "+s);
				traceFileName = s+".txt";
				
				// read fault trace from file
				file = new File(ROOT_DATA_DIR+FAULT_TRACE_DIR_NAME+traceFileName);
				for (String line : Files.readLines(file, Charset.defaultCharset())) {
					//			System.out.println(line);
					line = line.trim();
					String[] split = line.split("\t");	// tab delimited
					Preconditions.checkState(split.length == 2, "Expected 2 items, got %s", split.length);
					//			System.out.println(split[0]+"\t"+split[1]);
					double lon = Double.valueOf(split[0]);
					double lat = Double.valueOf(split[1]);
					fltTrace.add(new Location(lat,lon,UPPER_SEIS_DEPTH));
				}

				FaultSectionPrefData fltSectData = new FaultSectionPrefData();
				fltSectData.setSectionId(s);
				fltSectData.setSectionName("Section "+s);
				fltSectData.setShortName("Sect"+s);
				fltSectData.setFaultTrace(fltTrace);
				fltSectData.setAveDip(FAULT_DIP);
				fltSectData.setAveUpperDepth(UPPER_SEIS_DEPTH);
				fltSectData.setAveLowerDepth(LOWER_SEIS_DEPTH);
				fltSectData.setAveSlipRate(sectSlipRate[s]);
				fltSectData.setSlipRateStdDev(sectSlipRateStdDev[s]);
				fltSectData.setAveRake(FAULT_RAKE);;
				
				if(D) {
					String str = new String();
					str += "sectionId = "+fltSectData.getSectionId()+"\n";
					str += "sectionName = "+fltSectData.getSectionName()+"\n";
					str += "shortName = "+fltSectData.getShortName()+"\n";
					str += "aveLongTermSlipRate = "+fltSectData.getOrigAveSlipRate()+"\n";
					str += "slipRateStdDev = "+fltSectData.getOrigSlipRateStdDev()+"\n";
					str += "aveDip = "+fltSectData.getAveDip()+"\n";
					str += "aveRake = "+fltSectData.getAveRake()+"\n";
					str += "aveUpperDepth = "+fltSectData.getOrigAveUpperDepth()+"\n";
					str += "aveLowerDepth = "+fltSectData.getAveLowerDepth()+"\n";
					str += "aseismicSlipFactor = "+fltSectData.getAseismicSlipFactor()+"\n";
					str += "couplingCoeff = "+fltSectData.getCouplingCoeff()+"\n";
					str += "dipDirection = "+fltSectData.getDipDirection()+"\n";
					str += "dateOfLastEventMillis = "+fltSectData.getDateOfLastEvent()+"\n";
					str += "slipInLastEvent = "+fltSectData.getSlipInLastEvent()+"\n";
					str += "traceLength = "+fltSectData.getTraceLength()+"\n";
					str += "downDipWidth = "+fltSectData.getOrigDownDipWidth()+"\n";
					str += "area (stirling surface) = "+fltSectData.getStirlingGriddedSurface(1.0).getArea()+"\n";
					str += "faultTrace:\n";
					for(int i=0; i <fltSectData.getFaultTrace().size();i++) {
						Location loc = fltSectData.getFaultTrace().get(i);
						str += "\t"+loc.getLatitude()+", "+loc.getLongitude()+", "+loc.getDepth()+"\n";
					}
					System.out.println(str);
				}
				
				faultSectionDataList.add(fltSectData);
				
			}
			

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public ArrayList<FaultSectionPrefData> getFaultSectionDataList() {
		return faultSectionDataList;
	}
	
	public ArrayList<SegRateConstraint> getSectionRateConstraints() {
		return sectionRateConstraints;
	}
	
	public int[][] getRupSectionMatrix() {
		return rupSectionMatrix;
	}
	
	/**
	 * @param args
	 */
	public static void main(String []args) {
		
		WasatchInversion wasatchInversion = new WasatchInversion();
		ArrayList<FaultSectionPrefData> fltSectDataList = wasatchInversion.getFaultSectionDataList();
		ArrayList<SegRateConstraint> sectionRateConstraints = wasatchInversion.getSectionRateConstraints();
		int[][] rupSectionMatrix = wasatchInversion.getRupSectionMatrix();

		//Write section data (e.g., so it can be checked)
		String str = "sectID\tslipRate\tslipRateStd\tlength\tdip\trake\tupDepth\tlowDepth\tDDW+\n";
		for(int s=0;s<fltSectDataList.size();s++) {
			FaultSectionPrefData fltSectData = fltSectDataList.get(s);
			str += fltSectData.getSectionId()+"\t";
			str += fltSectData.getOrigAveSlipRate()+"\t";
			str += fltSectData.getOrigSlipRateStdDev()+"\t";
			str += (float)fltSectData.getTraceLength()+"\t";
			str += fltSectData.getAveDip()+"\t";
			str += fltSectData.getAveRake()+"\t";
			str += fltSectData.getOrigAveUpperDepth()+"\t";
			str += fltSectData.getAveLowerDepth()+"\t";
			str += fltSectData.getOrigDownDipWidth()+"\t";
			str += (float)fltSectData.getStirlingGriddedSurface(1.0).getArea()+"\n";
		}
		System.out.println(str);

//		System.exit(0);
		
		// this will be used to keep track of runtimes
		long startTimeMillis = System.currentTimeMillis();
		if(D)
			System.out.println("Starting Inversion");

		// set inversion attributes
//		String slipModelType = FaultSystemRuptureRateInversion.UNIFORM_SLIP_MODEL;
		String slipModelType = FaultSystemRuptureRateInversion.TAPERED_SLIP_MODEL;
		MagAreaRelationship magAreaRel = new HanksBakun2002_MagAreaRel();
//		MagAreaRelationship magAreaRel = new Ellsworth_B_WG02_MagAreaRel();
		double relativeSectRateWt=1;
		double relative_aPrioriRupWt = 0;	// 
		String aPrioriRupRatesFilename = ROOT_PATH+"data/aPrioriRupRatesFromGR_MFD.txt";
		double relative_smoothnessWt = 0;	// KEEP ZERO UNTIL THIS IS PROPERLY IMPLEMENTED
		boolean wtedInversion = true;
		double minRupRate = 1e-8;
		boolean applyProbVisible = true;
		double moRateReduction = 0.1;	// this is the amount to reduce due to smaller earthquakes being ignored (not due to asiesmity or coupling coeff, which are part of the fault section attributes)
		double relativeMFD_constraintWt = 0; // setting this to 1e6
		
		// GR Constraint (values obtained from first running FaultSystemRuptureRateInversion.getGR_DistFit() with relativeMFD_constraintWt=0)
		double minMag=6.4;
		double maxMag=8.1;
		double delta = 0.1;
		int num = 18;
		double moRate=2.415E17;
		GutenbergRichterMagFreqDist grConstraint = new GutenbergRichterMagFreqDist(minMag,num,delta);
		grConstraint.setAllButTotCumRate(minMag, maxMag, moRate, 1.0);

		// run the inversion
		FaultSystemRuptureRateInversion fltSysRupInversion = new  FaultSystemRuptureRateInversion(
				fltSectDataList, 
				sectionRateConstraints, 
				rupSectionMatrix, 
				slipModelType, 
				magAreaRel, 
				relativeSectRateWt, 
				relative_aPrioriRupWt, 
				aPrioriRupRatesFilename,
				relative_smoothnessWt, 
				wtedInversion, 
				minRupRate, 
				applyProbVisible, 
				moRateReduction,
				grConstraint,
				relativeMFD_constraintWt);
		
		// Non-negative least squares
//		fltSysRupInversion.doInversionNNLS();
		
		// Simulated annealing
		long numIterations = (long) 1e5;
		boolean initStateFromAprioriRupRates = false;
		fltSysRupInversion.doInversionSA(numIterations, initStateFromAprioriRupRates);


		double runTimeSec = ((double)(System.currentTimeMillis()-startTimeMillis))/1000.0;
		System.out.println("Done with Inversion after "+(float)runTimeSec+" seconds.");
		
		// do this once to write out a priori rupture rates from MFD constraint
//		fltSysRupInversion.writeApriorRupRatesFromMFD_Constrint(ROOT_PATH+"data/aPrioriRupRatesFromGR_MFD.txt");
		
		// these write to system out, not to a file
		fltSysRupInversion.writeFinalStuff();
		fltSysRupInversion.writePredErrorInfo();

		String dirName = ROOT_PATH+"OutputFigsAndData";
	    File file = new File(dirName);
	    file.mkdirs();
	    fltSysRupInversion.plotStuff(dirName);
	    fltSysRupInversion.plotMagHistograms();
	    fltSysRupInversion.writeAndPlotSegPartMFDs(dirName, true);
	    fltSysRupInversion.writeAndPlotNonZeroRateRups(dirName, true);
		
	}



}
