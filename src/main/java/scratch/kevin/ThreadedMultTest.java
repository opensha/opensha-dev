package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.StatUtils;

import scratch.UCERF3.utils.MatrixIO;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;

public class ThreadedMultTest extends Thread {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/test_inversion/bench/agu/ncal_constrained");
		DoubleMatrix2D A = MatrixIO.loadSparse(new File(dir, "A.bin"), SparseCCDoubleMatrix2D.class);
		DoubleMatrix1D d = new DenseDoubleMatrix1D(MatrixIO.doubleArrayFromFile(new File(dir, "d.bin")));
		
		for (int i=32; i<34; i++)
			doTest(i, A, d);
	}
	
	private static void doTest(int num, DoubleMatrix2D A, DoubleMatrix1D d) throws InterruptedException {
		System.out.println("*** "+num+" THREADS ***");
		
		ArrayList<ThreadedMultTest> threads = new ArrayList<ThreadedMultTest>();
		
		for (int i=0; i<num; i++)
			threads.add(new ThreadedMultTest(A, d));
		
		for (Thread t : threads)
			t.start();
		
		for (Thread t : threads)
			t.join();
		
		double[] mults = new double[num];
		long[] times = new long[num];
		double[] mps = new double[num];
		long maxTime = 0;
		
		for (int i=0; i<num; i++) {
			ThreadedMultTest t = threads.get(i);
			mults[i] = t.mults;
			times[i] = t.totTime;
			if (t.totTime > maxTime)
				maxTime = t.totTime;
			mps[i] = (double)t.mults / (t.totTime / 1000d);
		}
		
		double overallMPS = StatUtils.sum(mults) / (maxTime / 1000d);
		
		System.out.println("Total mults: "+StatUtils.sum(mults));
		System.out.println("AVG mults: "+StatUtils.mean(mults));
		System.out.println("AVG mults/sec: "+StatUtils.mean(mps));
		System.out.println("OVERALL mults/sec: "+overallMPS);
		System.out.println("*******************************");
		System.out.println();
	}
	
	private DoubleMatrix2D A;
	private DoubleMatrix1D d;
	private DoubleMatrix1D x;
	
	private int mults = 0;
	private long totTime;
	
	public ThreadedMultTest(DoubleMatrix2D A, DoubleMatrix1D d) {
		this.A = A;
		this.d = d;
		x = new DenseDoubleMatrix1D(A.columns());
	}

	@Override
	public void run() {
		StopWatch watch = new StopWatch();
		
		watch.start();
		while (watch.getTime() < 15000) {
			A.zMult(x, d);
			mults++;
		}
		watch.stop();
		totTime = watch.getTime();
	}

}
