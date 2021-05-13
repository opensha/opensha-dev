package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

import scratch.UCERF3.utils.MatrixIO;

public class SparseWriteTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("D:\\Documents\\temp\\UCERF3_Geologic_Model1");
		File input = new File(dir, "A.bin");
		File output = new File(dir, "output.bin");
//		input = output;
		System.out.println("Loading...");
		Stopwatch watch = Stopwatch.createStarted();
//		DoubleMatrix2D mat = MatrixIO.loadSparse(input, SparseCCDoubleMatrix2D.class);
		DoubleMatrix2D mat = MatrixIO.loadSparse(input, SparseDoubleMatrix2D.class);
		watch.stop();
		System.out.println("Loading took "+watch.elapsed(TimeUnit.SECONDS)+" secs");
		System.out.println("Saving...");
		watch.reset();
		watch.start();
		MatrixIO.saveSparse(mat, output);
		watch.stop();
		System.out.println("Saving took "+watch.elapsed(TimeUnit.SECONDS)+" secs");
		System.out.println("DONE.");
	}

}
