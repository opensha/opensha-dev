package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import static org.junit.Assert.*;

import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.junit.Test;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.aftershockStatisticsETAS.griddedInterpGMPE.NDimArrayCalc;
import scratch.aftershockStatisticsETAS.griddedInterpGMPE.NDimensionalLinearInterpolation;

public class NDimensionalLinearInterpolationTest {
	
	private static Random r = new Random();
	
	private static final double tolerance = 1e-5;;
	
	/*
	 * test index retrieval, no actual interpolation
	 */
	
	@Test
	public void testRetrieval1D() {
		validateRetrieval(new int[] { 5+r.nextInt(5) });
	}
	
	@Test
	public void testRetrieval2D() {
		validateRetrieval(new int[] { 5+r.nextInt(5), 5+r.nextInt(5) });
	}
	
	@Test
	public void testRetrieval3D() {
		validateRetrieval(new int[] { 2+r.nextInt(3), 2+r.nextInt(3), 2+r.nextInt(3) });
	}
	
	@Test
	public void testRetrieval4D() {
		validateRetrieval(new int[] { 2+r.nextInt(1), 2+r.nextInt(1), 2, 2 });
	}
	
	private void validateRetrieval(int[] dimensions) {
		NDimArrayCalc arrayCalc = new NDimArrayCalc(dimensions);
		double[] flatArray = new double[arrayCalc.rawArraySize()];
		for (int i=0; i<flatArray.length; i++)
			flatArray[i] = r.nextDouble();
		validateRetrieval(flatArray, arrayCalc);
	}
	
	private void validateRetrieval(double[] flatArray, NDimArrayCalc arrayCalc) {
		List<int[]> indexesList = NDimensionalLinearInterpolation.getSubInterpIndexes(arrayCalc.getDimensions());
		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation(arrayCalc.getNumDimensions());
		System.out.println("Validating retrieval for "+indexesList.size()+" indexes");
		for (int[] indexes : indexesList) {
			int flatIndex = arrayCalc.getIndex(indexes);
			// as a double
			double[] interpIndexes = new double[indexes.length];
			for (int i=0; i<indexes.length; i++)
				interpIndexes[i] = indexes[i];
			double expected = flatArray[flatIndex];
			double actual = interp.interpolate(flatArray, arrayCalc, interpIndexes);
			assertEquals(expected, actual, tolerance);
		}
	}
	
	
	/*
	 * test actual interpolation
	 */

	@Test
	public void test1D_forward() {
		int numTests = 10;
		int numEach = 10;
		doTest1D(numTests, numEach, false);
	}

	@Test
	public void test1D_reversed() {
		int numTests = 10;
		int numEach = 10;
		doTest1D(numTests, numEach, true);
	}
	
	private void doTest1D(int numTests, int numEach, boolean reversed) {
		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation(1);
		for (int t=0; t<numTests; t++) {
			double min = r.nextDouble() * 1000 - 500;
			double max = min + r.nextDouble() * 7500;
			double[] data;
			if (reversed)
				data = new double[] { max, min };
			else
				data = new double[] { min, max };
			for (int i=0; i<numEach; i++) {
				double fract = (double)i/(double)(numEach-1);
				double expected = ref1D(data, fract);
				double[] deltas = { fract };
				double interpVal = interp.interpolateTrimmed(data, deltas);
				assertEquals("i="+i+", fract="+fract+", data=["+data[0]+","+data[1]+"]", expected, interpVal, tolerance);
			}
		}
	}
	
	private static double ref1D(double[] data, double x) {
		int x0 = (int)x;
		if (x0 == x)
			return data[x0];
		int x1 = x0+1;
		Preconditions.checkState(x0 >= 0 && x1 < data.length, "Bad index. Len=%s, x=%s, x0=%s, x1=%s", data.length, x, x0, x1);
		double delta = x - x0;
		return (1d - delta)*data[x0] + delta*data[x1];
	}
	
	@Test
	public void test2D_simple_2x2() {
		double[][] data = new double[2][2];
		data[0][0] = 0d;
		data[1][1] = 1d;
		data[0][1] = 0.5;
		data[1][0] = -0.5;
		
		doTest2D(data);
	}
	
	@Test
	public void test2D_NxN() {
		int numTests = 20;
		for (int t=0; t<numTests; t++) {
			double[][] data = new double[2+r.nextInt(30)][2+r.nextInt(30)];
			for (int x=0; x<data.length; x++)
				for (int y=0; y<data[x].length; y++)
					data[x][y] = r.nextDouble();
			
			doTest2D(data);
		}
	}
	
	private void doTest2D(double[][] data) {
		NDimArrayCalc arrayCalc = new NDimArrayCalc(new int[] { data.length, data[0].length });
		double[] flatArray = new double[arrayCalc.rawArraySize()];
		
		for (int x=0; x<data.length; x++)
			for (int y=0; y<data[x].length; y++)
				flatArray[arrayCalc.getIndex(x, y)] = data[x][y];
		
		int numTestPerDimension = 20;
		
		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation(2);
		
		for (int i=0; i<numTestPerDimension; i++) {
			double x = ((double)i/(double)(numTestPerDimension-1)) * (data.length-1);
			for (int j=0; j<numTestPerDimension; j++) {
				double y = ((double)j/(double)(numTestPerDimension-1)) * (data[0].length-1);
				double expected = ref2D(data, x, y);
				double actual = interp.interpolate(flatArray, arrayCalc, new double[] { x, y });
				assertEquals("x["+i+"]="+x+", y["+j+"]="+y, expected, actual, tolerance);
			}
		}
	}
	
	@Test
	public void test3D_simple_2x2x2() {
		double[][][] data = new double[2][2][2];
		data[0][0][0] = 0d;
		data[1][1][0] = 1d;
		data[0][1][0] = 0.5;
		data[1][0][0] = -0.5;
		data[0][0][1] = 1d;
		data[1][1][1] = 2d;
		data[0][1][1] = 1.5;
		data[1][0][1] = 0.5;
		
		doTest3D(data);
	}
	
	@Test
	public void test3D_NxNxN() {
		int numTests = 20;
		for (int t=0; t<numTests; t++) {
			double[][][] data = new double[2+r.nextInt(30)][2+r.nextInt(30)][2+r.nextInt(30)];
			for (int x=0; x<data.length; x++)
				for (int y=0; y<data[x].length; y++)
					for (int z=0; z<data[x][y].length; z++)
						data[x][y][z] = r.nextDouble();
			
			doTest3D(data);
		}
	}
	
	private void doTest3D(double[][][] data) {
		NDimArrayCalc arrayCalc = new NDimArrayCalc(new int[] { data.length, data[0].length, data[0][0].length });
		double[] flatArray = new double[arrayCalc.rawArraySize()];
		
		for (int x=0; x<data.length; x++)
			for (int y=0; y<data[x].length; y++)
				for (int z=0; z<data[x][y].length; z++)
					flatArray[arrayCalc.getIndex(x, y, z)] = data[x][y][z];
		
		int numTestPerDimension = 20;
		
		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation(3);
		
		for (int i=0; i<numTestPerDimension; i++) {
			double x = ((double)i/(double)(numTestPerDimension-1)) * (data.length-1);
			for (int j=0; j<numTestPerDimension; j++) {
				double y = ((double)j/(double)(numTestPerDimension-1)) * (data[0].length-1);
				for (int k=0; k<numTestPerDimension; k++) {
					double z = ((double)k/(double)(numTestPerDimension-1)) * (data[0][0].length-1);
					double expected = ref3D(data, x, y, z);
					double actual = interp.interpolate(flatArray, arrayCalc, new double[] { x, y, z });
					assertEquals("x="+x+", y="+y+", z="+z, expected, actual, tolerance);
				}
			}
		}
	}
	
	private static double ref2D(double[][] data, double x, double y) {
		// from https://en.wikipedia.org/wiki/Bilinear_interpolation#Nonlinear
		int x0 = (int)x;
		int y0 = (int)y;
		int x1 = x0+1;
		int y1 = y0+1;
		
		double dx = x - x0;
		double dy = y - y0;
		
		if (x == x0)
			return ref1D(data[x0], y);
		else if (y == y0)
			return ref1D(new double[] {data[x0][y0], data[x1][y0]}, dx);
		
		double a00 = data[x0][y0];
		double a10 = data[x1][y0] - data[x0][y0];
		double a01 = data[x0][y1] - data[x0][y0];
		double a11 = data[x1][y1] + data[x0][y0] - (data[x1][y0] + data[x0][y1]);
		return a00  + a10*dx + a01*dy + a11*dx*dy;
	}
	
	private static double ref3D(double[][][] data, double x, double y, double z) {
		// from https://en.wikipedia.org/wiki/Bilinear_interpolation#Nonlinear
		int x0 = (int)x;
		int y0 = (int)y;
		int z0 = (int)z;
		int z1 = z0+1;
		
		double dz = z - z0;
		
		if (x == x0)
			return ref2D(data[x0], y, z);
		else if (y == y0) {
			double[][] subData = new double[data.length][data[0][0].length];
			for (int i=0; i<subData.length; i++)
				for (int j=0; j<subData[i].length; j++)
					subData[i][j] = data[i][y0][j];
			return ref2D(subData, x, z);
		} else if (z == z0) {
			double[][] subData = new double[data.length][data[0].length];
			for (int i=0; i<subData.length; i++)
				for (int j=0; j<subData[i].length; j++)
					subData[i][j] = data[i][j][z0];
			return ref2D(subData, x, y);
		}
		
		double[][] dataZ0 = new double[data.length][data[0].length];
		double[][] dataZ1 = new double[data.length][data[0].length];
		
		for (int i=0; i<data.length; i++) {
			for (int j=0; j<data[0].length; j++) {
				dataZ0[i][j] = data[i][j][z0];
				dataZ1[i][j] = data[i][j][z1];
			}
		}
		
		double interpZ0 = ref2D(dataZ0, x, y);
		double interpZ1 = ref2D(dataZ1, x, y);
		
		return ref1D(new double[] { interpZ0, interpZ1 }, dz);
	}
	
	@Test
	public void test3Dspeed() {
		int nx = 20;
		int ny = 36;
		int nz = 20;
		double[][][] data = new double[nx][ny][nx];
//		NDimArrayCalc arrayCalc = new NDimArrayCalc(new int[] { nx, ny, nz});
		NDimArrayCalc arrayCalc = new NDimArrayCalc(new int[] { nx, ny, nz, 2 });
		double[] flatData = new double[arrayCalc.rawArraySize()];
		for (int x=0; x<nx; x++) {
			for (int y=0; y<ny; y++) {
				for (int z=0; z<nz; z++) {
					data[x][y][z] = r.nextDouble();
//					flatData[arrayCalc.getIndex(x, y, z)] = data[x][y][z];
					flatData[arrayCalc.getIndex(x, y, z, 0)] = data[x][y][z];
					flatData[arrayCalc.getIndex(x, y, z, 1)] = data[x][y][z];
				}
			}
		}
		
		NDimensionalLinearInterpolation interp = new NDimensionalLinearInterpolation(arrayCalc.getNumDimensions());
		
		int numTests = 1000000;
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<numTests; i++) {
			double x = r.nextDouble()*(nx-1);
			double y = r.nextDouble()*(ny-1);
			double z = r.nextDouble()*(nz-1);
			ref3D(data, x, y, z);
//			interp.interpolate(flatData, arrayCalc, new double[] {x,y,z});
		}
		watch.stop();
		double secs = (double)watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Ref took "+secs+" secs");
		
		watch = Stopwatch.createStarted();
		for (int i=0; i<numTests; i++) {
			double x = r.nextDouble()*(nx-1);
			double y = r.nextDouble()*(ny-1);
			double z = r.nextDouble()*(nz-1);
//			ref3D(data, x, y, z);
//			interp.interpolate(flatData, arrayCalc, new double[] {x,y,z});
			interp.interpolate(flatData, arrayCalc, new double[] {x,y,z,0});
		}
		watch.stop();
		secs = (double)watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("NDim took "+secs+" secs");
	}

}
