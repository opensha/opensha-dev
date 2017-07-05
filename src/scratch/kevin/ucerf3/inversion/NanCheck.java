package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;

import com.google.common.base.Preconditions;

import scratch.UCERF3.utils.MatrixIO;

public class NanCheck {
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/tmp/FM3_1_GEOB_MaAvU2_DsrTap_DrAvU2_Char_inputs");
		File a = new File(dir, "a.bin");
		MatrixIO.loadSparse(a);
		File a_ineq = new File(dir, "a_ineq.bin");
		MatrixIO.loadSparse(a_ineq);
		File d = new File(dir, "d.bin");
		for (double val : MatrixIO.doubleArrayFromFile(d))
			Preconditions.checkState(!Double.isNaN(val));
		File d_ineq = new File(dir, "d_ineq.bin");
		for (double val : MatrixIO.doubleArrayFromFile(d_ineq))
			Preconditions.checkState(!Double.isNaN(val));
		File initial = new File(dir, "initial.bin");
		for (double val : MatrixIO.doubleArrayFromFile(initial))
			Preconditions.checkState(!Double.isNaN(val));
	}

}
