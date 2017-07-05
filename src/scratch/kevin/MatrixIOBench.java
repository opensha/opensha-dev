package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.util.ClassUtils;

import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Stopwatch;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class MatrixIOBench {
	
	private static float getSecs(Stopwatch watch) {
		return (float)(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d);
	}
	
	private static void ioBench(File a_file, Class<? extends DoubleMatrix2D> clazz) throws IOException {
		String cname = ClassUtils.getClassNameWithoutPackage(clazz);
		System.out.println(cname+": loading");
		Stopwatch watch = Stopwatch.createStarted();
		DoubleMatrix2D mat = MatrixIO.loadSparse(a_file, clazz);
		watch.stop();
		System.out.println(cname+": took "+getSecs(watch)+" secs to load.");
		
		watch.reset();
		System.out.println(cname+": saving");
		File tempFile = File.createTempFile("openSHA", "a_matrix.bin");
		watch.start();
		MatrixIO.saveSparse(mat, tempFile);
		watch.stop();
		System.out.println(cname+": took "+getSecs(watch)+" secs to save.");
		
		tempFile.delete();
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File a_file = new File("/home/kevin/OpenSHA/UCERF3/inversions/2012_02_22-UCERF3_Geologic_Model2/a.bin");
		
		System.out.println("Start your profiling!");
		try {
			Thread.sleep(5000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ioBench(a_file, SparseDoubleMatrix2D.class);
		System.gc();
		ioBench(a_file, SparseCCDoubleMatrix2D.class);
	}

}
