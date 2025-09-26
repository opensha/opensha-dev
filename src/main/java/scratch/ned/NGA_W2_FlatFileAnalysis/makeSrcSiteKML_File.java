package scratch.ned.NGA_W2_FlatFileAnalysis;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;

import scratch.UCERF3.analysis.FaultBasedMapGen;

public class makeSrcSiteKML_File {

	public static void main(String[] args) {
		
        // read csv rate file
		String csvPath = "/Users/field/Library/CloudStorage/OneDrive-DOI/Field_Other/MiscDocs/test.csv";
		double minMagThresh=4;
		double maxMagThresh=9;
		double maxDist = 100;
				
		File saveDir = new File("TEST_DIR_HERE");
		if(!saveDir.exists()) saveDir.mkdir();

		File file = new File(csvPath);
		CSVFile<String> csvFile=null;
		try {
			csvFile = CSVFile.readFile(file, true);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
//		CPT cpt = new CPT(minMagThresh, maxMagThresh, Color.BLUE, Color.RED);
		CPT cpt = new CPT(4.0, 9.0, Color.BLACK, Color.BLUE, Color.GREEN, Color.ORANGE, Color.RED, Color.MAGENTA);
//		GMT_CPT_Files cptFile = GMT_CPT_Files.MAX_SPECTRUM;
//		CPT cpt=null;
//		try {
//			cpt = cptFile.instance();
//		} catch (IOException e1) {
//			e1.printStackTrace();
//		}
//		cpt.rescale(4, 9);


			
			// sort file rows by mag
//			ArrayList<Integer> csvRowIndicesSortedByMag = new ArrayList<Integer>();
//			double deltaMag = 0.25;
//			for(double curMag=minMagThresh; curMag<maxMagThresh; curMag+=deltaMag) {
//				for(int i=1; i<csvFile.getNumRows(); i++) {
//					double mag = csvFile.getDouble(i, 2);
//					if(mag>=curMag && mag<curMag+deltaMag)
//						csvRowIndicesSortedByMag.add(i);
//				}
//			}
////			Collections.reverse(csvRowIndicesSortedByMag);
//			System.out.println(csvFile.getDouble(csvRowIndicesSortedByMag.get(0), 2));
//			System.out.println(csvFile.getDouble(csvRowIndicesSortedByMag.get(csvRowIndicesSortedByMag.size()-1), 2));
			
		for(double magMid=5; magMid<8.25; magMid+=1.0) {
			double minMag = magMid-0.5;
			double maxMag = magMid+0.5;
			
			List<LocationList> pathList = new ArrayList<LocationList>();
			List<String> pathNameList = new ArrayList<String>();
			List<Double> magList = new ArrayList<Double>();
			int numPaths=0;


			
			for(int i=1; i<csvFile.getNumRows(); i++) {
				double mag = csvFile.getDouble(i, 2);
				if(mag<minMag || mag>=maxMag)
					continue;
				String qkName = csvFile.get(i, 0);
				String stationName = csvFile.get(i, 1);
				if(mag >= minMagThresh && mag <= maxMagThresh) {
					double hypoLat = csvFile.getDouble(i, 3);
					double hypoLon = csvFile.getDouble(i, 4);
					double hypoDepth = csvFile.getDouble(i, 5);
					double stationLat = csvFile.getDouble(i, 6);
					double stationLon = csvFile.getDouble(i, 7);
					String pathName = qkName+"-->"+stationName;

					if(hypoLat>90 || hypoLat<-90) {
						System.out.println("SKIPPING: hypoLat="+hypoLat+"\tfor "+pathName);
						continue;
					}
					if(stationLat>90 || stationLat<-90) {
						System.out.println("SKIPPING: stationLat="+stationLat+"\tfor "+pathName);
						continue;
					}

					if(hypoLon>360 || hypoLon<-360) {
						System.out.println("SKIPPING: hypoLon="+hypoLon+"\tfor "+pathName);
						continue;
					}
					if(stationLon>360 || stationLon<-360) {
						System.out.println("SKIPPING: stationLon="+stationLon+"\tfor "+pathName);
						continue;
					}

					Location hypoLoc = new Location(hypoLat,hypoLon);
					Location stationLoc = new Location(stationLat,stationLon);

					if(LocationUtils.horzDistance(hypoLoc, stationLoc)> maxDist)
						continue;;


						LocationList locList = new LocationList();
						locList.add(hypoLoc);
						locList.add(stationLoc);


						pathList.add(locList);
						pathNameList.add(pathName);
						magList.add(mag);

						numPaths+=1;
						//					System.out.println(pathName+"\t"+mag);
				}
			}
			System.out.println("numPaths="+numPaths+" for magMid="+magMid);

			double[] magArray = new double[magList.size()];
			for(int i=0;i<magList.size();i++)
				magArray[i] = magList.get(i);
			
			
			String prefix = "testFileHere_"+magMid;
			boolean skipNans = true;
			int numColorBins=20;
			int lineWidth=2;
			String name = "magMid_"+magMid;
			
			// scratch.UCERF3.analysis.FaultBasedMapGen.makeFaultKML(CPT, List<LocationList>, double[], File, String, boolean, int, int, String, List<String>)
			try {
				FaultBasedMapGen.makeFaultKML(cpt, pathList, magArray,saveDir, prefix, skipNans, numColorBins, 
						lineWidth, name, pathNameList);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		
		}
		
			
		
//		public static void makeFaultKML(CPT cpt, List<LocationList> faults, double[] values,
//				File saveDir, String prefix, boolean skipNans, int numColorBins, int lineWidth,
//				String name, List<String> descriptions)

	}

}
