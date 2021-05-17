package scratch.kevin;

import java.io.FileWriter;
import java.io.IOException;

import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.ArbDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.XYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.ClassUtils;

public class GeoDataSetBenchmark {
	
	private static double[] benchmark(XYZ_DataSet xyz, int num) {
		long putStart = System.currentTimeMillis();
		double inc = 0.0001;
		double lat = 0;
		double lon = 0;
		for (int i=0; i<num; i++) {
			xyz.set(lat, lon, i);
			lat += inc;
			lon += inc;
		}
		long putEnd = System.currentTimeMillis();

		lat = 0;
		lon = 0;
		long getStart = System.currentTimeMillis();
		for (int i=0; i<num; i++) {			
			xyz.get(lat, lon);
			lat += inc;
			lon += inc;
		}
		long getEnd = System.currentTimeMillis();
		
		double setSecs = ((double)(putEnd - putStart))/1000d;
		double getSecs = ((double)(getEnd - getStart))/1000d;
		
		System.out.println("******* "+ClassUtils.getClassNameWithoutPackage(xyz.getClass()));
		System.out.println("set: "+(float)setSecs + " s");
		System.out.println("get: "+(float)getSecs + " s");
		
		double[] ret = { setSecs, getSecs };
		return ret;
	}
	
	private static void sleepGC() {
		System.gc();
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
}
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		FileWriter fw = new FileWriter("/home/kevin/OpenSHA/xyz_bench/compare.csv");
//		fw.write(",SET,,,,GET,,,\n");
//		fw.write("size,ArbDiscrXYZ_DataSet,ArbDiscrGeoDataSet,ArrayListBasedGeoDataSet,ArrayBasedGeoDataSet," +
//				"ArbDiscrXYZ_DataSet,ArbDiscrGeoDataSet,ArrayListBasedGeoDataSet,ArrayBasedGeoDataSet\n");
		int adder = 20000;
		for (int size=20000; size<=220000; size+=adder) {
			sleepGC();
			System.out.println("SIZE: "+size);
//			double[] results0 = benchmark(new ArbDiscrXYZ_DataSet(), size);
//			sleepGC();
			double[] results1 = benchmark(new ArbDiscrGeoDataSet(true), size);
			sleepGC();
//			double[] results2 = benchmark(new ArrayListBasedGeoDataSet(true), size);
//			sleepGC();
//			double[] results3 = benchmark(new ArrayBasedGeoDataSet(true, size), size);
//			sleepGC();
//			
//			fw.write(size+","+results0[0]+","+results1[0]+","+results2[0]+","+results3[0]+","+
//					results0[1]+","+results1[1]+","+results2[1]+","+results3[1]+"\n");
//			fw.flush();
			if (size == 100000)
				adder = 40000;
		}
//		fw.close();
	}

}

