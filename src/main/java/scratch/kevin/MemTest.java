package scratch.kevin;

import java.text.DecimalFormat;
import java.util.Random;

import com.google.common.base.Preconditions;

public class MemTest {

	public static void main(String[] args) {
		long maxMemBytes = Runtime.getRuntime().maxMemory();
		long gbToBytes = 1024l*1024l*1024l;
		double maxMemGB = (double)maxMemBytes/(double)gbToBytes;
		long curMemBytes = Runtime.getRuntime().totalMemory();
		double curMemGB = (double)curMemBytes/(double)gbToBytes;
		DecimalFormat gbDF = new DecimalFormat("0.00");
		System.gc();
		System.out.println("Max VM memory: "+gbDF.format(maxMemGB)+" GB");
		System.out.println("Starting VM memory: "+gbDF.format(curMemGB)+" GB");
		int targetGB = Integer.parseInt(args[0]);
		long targetBytes = gbToBytes*targetGB;
		long numLongs = targetBytes / 8l;
		System.out.println("Will try to allocate "+targetGB+" GB = "+targetBytes+" bytes = "+numLongs+" longs");
		Preconditions.checkState(numLongs < (long)Integer.MAX_VALUE);
		long[] hugeArray = new long[(int)numLongs];
		System.out.println("Done initializing array");
		curMemBytes = Runtime.getRuntime().totalMemory();
		curMemGB = (double)curMemBytes/(double)gbToBytes;
		System.out.println("Current VM memory: "+gbDF.format(curMemGB));
		System.gc();
		Random r = new Random();
		while (true)
			hugeArray[r.nextInt(hugeArray.length)] = r.nextLong();
	}

}
