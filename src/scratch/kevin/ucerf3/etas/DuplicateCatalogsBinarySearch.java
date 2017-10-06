package scratch.kevin.ucerf3.etas;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;

public class DuplicateCatalogsBinarySearch {
	
	public static final int buffer_len = 655360;
	public static final int bytes_per_rup = 70;

	public static void main(String[] args) throws IOException, NoSuchAlgorithmException {
		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m5_preserve.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_22-mojave_m7-10yr-full_td-no_ert-combined/results_m5_preserve.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_12_03-mojave_m7-10yr-gridded-only/results_descendents_m5_preserve.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2017_01_02-haywired_m7-10yr-gridded-only-200kcombined/results_descendents_m5_preserve.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_06_15-haywired_m7-10yr-full_td-no_ert-combined/results_descendents_m5.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_08_30-san_jacinto_0_m4p8-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined/results_descendents.bin");
//		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2017_04_13-parkfield-10yr-full_td-no_ert-combined/results_descendents_m4_preserve.bin");
		
		FileInputStream fis = new FileInputStream(catalogFile);
		BufferedInputStream bis = new BufferedInputStream(fis, buffer_len);
		
		DataInputStream in = new DataInputStream(bis);
		
		int numCatalogs = in.readInt();
		
		int numEmpty = 0;
		int numDups = 0;
		Map<String, List<Integer>> values = new HashMap<>();
		Joiner j = Joiner.on(",");
		int maxIdentical = 0;
		for (int i=0; i<numCatalogs; i++) {
			MessageDigest md = MessageDigest.getInstance("MD5");
			in.readShort(); // version
			int numRups = in.readInt();
//			System.out.println(i+" Rups: "+numRups);
			Preconditions.checkState(numRups >= 0, "Bad numRups: %s", numRups);
			if (numRups == 0) {
				numEmpty++;
				continue;
			}
			int numBytesForCat = bytes_per_rup*numRups;
			byte[] catBytes = new byte[numBytesForCat];
			int read = in.read(catBytes);
			Preconditions.checkState(read == numBytesForCat, "Bad read. Expected %s, read %s", numBytesForCat, read);
			byte[] digest = md.digest(catBytes);
			String digestStr = bytesToHex(digest);
//			System.out.println(i+" Digest: "+digestStr);
			if (values.containsKey(digestStr)) {
				System.out.println(i+" is a duplicate of "+j.join(values.get(digestStr)));
				numDups++;
				values.get(digestStr).add(i);
				int num = values.get(digestStr).size();
				if (num > maxIdentical)
					maxIdentical = num;
			} else {
				List<Integer> list = new ArrayList<>();
				list.add(i);
				values.put(digestStr, list);
			}
		}
		
		System.out.println(numDups+"/"+numCatalogs+" duplicates (excluding empty)");
		System.out.println(numEmpty+"/"+numCatalogs+" empties");
		System.out.println("Max identical: "+maxIdentical);
		
		in.close();
	}
	
	private final static char[] hexArray = "0123456789ABCDEF".toCharArray();
	public static String bytesToHex(byte[] bytes) {
	    char[] hexChars = new char[bytes.length * 2];
	    for ( int j = 0; j < bytes.length; j++ ) {
	        int v = bytes[j] & 0xFF;
	        hexChars[j * 2] = hexArray[v >>> 4];
	        hexChars[j * 2 + 1] = hexArray[v & 0x0F];
	    }
	    return new String(hexChars);
	}

}
