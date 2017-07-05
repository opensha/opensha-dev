package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.util.DataUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.utils.MatrixIO;

public class RatesFileCompare {

	public static void main(String[] args) throws IOException {
		File file1 = new File("/tmp/asdf/rates.bin");
		File file2 = new File("/tmp/td_output_POISSON.bin");
		
		double[] rates1 = MatrixIO.doubleArrayFromFile(file1);
		double[] rates2 = MatrixIO.doubleArrayFromFile(file2);
		
		Preconditions.checkState(rates1.length == rates2.length);
		
		int num = rates1.length;
		
		int numIdentical = 0;
		int numIdenticalFloat = 0;
		int numClose = 0;
		int numZero1 = 0;
		int numZero2 = 0;
		
		for (int i=0; i<num; i++) {
			double v1 = rates1[i];
			double v2 = rates2[i];
			
			if (v1 == v2)
				numIdentical++;
			if ((float)v1 == (float)v2)
				numIdenticalFloat++;
			if (DataUtils.getPercentDiff(v2, v1) < 0.01)
				numClose++;
			if (v1 == 0d)
				numZero1++;
			if (v2 == 0d)
				numZero2++;
		}
		
		System.out.println(numIdentical+"/"+num+" identical double precision");
		System.out.println(numIdenticalFloat+"/"+num+" identical float precision");
		System.out.println(numClose+"/"+num+" within 0.01%");
		System.out.println(numZero1+"/"+num+" zero from file 1");
		System.out.println(numZero2+"/"+num+" zero from file 2");
	}

}
