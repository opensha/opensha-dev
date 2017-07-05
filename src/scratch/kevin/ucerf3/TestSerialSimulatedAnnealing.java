package scratch.kevin.ucerf3;

import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Before;
import org.junit.Test;

import scratch.UCERF3.simulatedAnnealing.SerialSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.ThreadedSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.completion.IterationCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.TimeCompletionCriteria;

import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class TestSerialSimulatedAnnealing {
	
	private SparseDoubleMatrix2D a_reg;
	private SparseDoubleMatrix2D a_ineq_reg;
	private SparseCCDoubleMatrix2D a_cc;
	private SparseCCDoubleMatrix2D a_ineq_cc;
	private double[] d;
	private double[] d_ineq;
	private double[] initialState;
	
	private double precision = 1e-10;
	
	private SerialSimulatedAnnealing sa_reg;
	private SerialSimulatedAnnealing sa_cc;

	@Before
	public void setUp() throws Exception {
		int rows = 1000;
		int cols = 10000;
		int numNonZero = 10000;
		a_reg = new SparseDoubleMatrix2D(rows, cols);
		a_cc = new SparseCCDoubleMatrix2D(rows, cols);
		
		Random r = new Random();
		for (int i=0; i<numNonZero; i++) {
			int row = r.nextInt(rows);
			int col = r.nextInt(cols);
			double val = r.nextDouble();
			
			a_reg.set(row, col, val);
			a_cc.set(row, col, val);
		}
		
		int ineq_rows = 100;
		int ineq_numNonZero = 1000;
		a_ineq_reg = new SparseDoubleMatrix2D(ineq_rows, cols);
		a_ineq_cc = new SparseCCDoubleMatrix2D(ineq_rows, cols);
		
		for (int i=0; i<ineq_numNonZero; i++) {
			int row = r.nextInt(ineq_rows);
			int col = r.nextInt(cols);
			double val = r.nextDouble();
			
			a_ineq_reg.set(row, col, val);
			a_ineq_cc.set(row, col, val);
		}

		d = new double[rows];
		for (int i=0; i<rows; i++)
			d[i] = Math.random();
		
		d_ineq = new double[ineq_rows];
		for (int i=0; i<ineq_rows; i++)
			d_ineq[i] = Math.random();
		
		initialState = new double[cols];
		for (int i=0; i<cols; i++)
			initialState[i] = Math.random();
		
		long seed = System.currentTimeMillis();
		Random r1 = new Random(seed);
		Random r2 = new Random(seed);
		
		sa_reg = new SerialSimulatedAnnealing(a_reg, d, initialState, 0, a_ineq_reg, d_ineq);
		sa_reg.setRandom(r1);
		sa_cc = new SerialSimulatedAnnealing(a_cc, d, initialState, 0, a_ineq_cc, d_ineq);
		sa_cc.setRandom(r2);
	}

	@Test
	public void testQuickMult() {
		long numIterations = 10000;
		
		System.out.println("Doing regular iterations...");
		sa_reg.iterate(numIterations);
		System.out.println("Doing column compressed iterations...");
		sa_cc.iterate(numIterations);
		System.out.println("Comparing results");
		
		System.out.println("Energy: "+sa_reg.getBestEnergy());
		
		verifySolutionsAndEnergy();
	}
	
	private void verifySolutionsAndEnergy() {
		assertEquals("energy doesn't match!", sa_reg.getBestEnergy()[0], sa_cc.getBestEnergy()[0], precision);
		
		double[] sol_reg = sa_reg.getBestSolution();
		double[] sol_cc = sa_cc.getBestSolution();
		
		for (int i=0; i<sol_reg.length; i++)
			assertEquals("solution doesn't match at index "+i, sol_reg[i], sol_cc[i], precision);
	}
	
	@Test
	public void testIndivVsBatch() {
		long numIterations = 100;
		
		long iter = 0l;
		
		for (int i=0; i<10000; i++) {
//			System.out.println("loop: "+i+" (energy: "+sa_cc.getBestEnergy()+")");
			sa_cc.iterate(iter, 0l, new IterationCompletionCriteria(iter+numIterations));
			iter += numIterations;
			
			double[] sol = sa_cc.getBestSolution();
			// now calculate energy the long way
			
			long seed = System.currentTimeMillis();
			sa_cc.setRandom(new Random(seed));
			sa_reg.setRandom(new Random(seed));
			
			IterationCompletionCriteria criteria = new IterationCompletionCriteria(iter + 1);
			
			// iterate each one for a single iteration
			// set energy to max value because one could keep it while the other doesn't due to double precision issues
			// in the energy compare step
			double[] initial = { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
					Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
			sa_reg.setResults(initial, sol);
			sa_cc.setResults(initial, sol);
			sa_reg.iterate(iter, 0l, criteria);
			sa_cc.iterate(iter, 0l, criteria);
			
			iter++;
			
			verifySolutionsAndEnergy();
		}
	}
	
	@Test
	public void testThreadedAcuracy() {
		ThreadedSimulatedAnnealing tsa = new ThreadedSimulatedAnnealing(a_cc, d, initialState, 0,
				a_ineq_cc, d_ineq, null, 4, TimeCompletionCriteria.getInSeconds(2));
		
		TimeCompletionCriteria criteria = TimeCompletionCriteria.getInSeconds(10);
		
		long iter = 0;
		for (int i=0; i<10; i++) {
			iter = tsa.iterate(iter, 0l, criteria)[0];
			// check that energy was computed correctly for the current solution
			
			double[] tsa_misfit = tsa.getBestMisfit();
			double[] tsa_ineq_misfit = tsa.getBestInequalityMisfit();
			
			SerialSimulatedAnnealing sa_test = new SerialSimulatedAnnealing(a_reg, d, tsa.getBestSolution(), 0d, a_ineq_reg, d_ineq);
			
			assertEquals("energy doesn't match!", sa_test.getBestEnergy()[0], tsa.getBestEnergy()[0], precision);
			
			double[] sa_misfit = sa_test.getBestMisfit();
			double[] sa_ineq_misfit = sa_test.getBestInequalityMisfit();
			
			for (int j=0; j<tsa_misfit.length; j++)
				assertEquals("solution doesn't match at index "+j, sa_misfit[j], tsa_misfit[j], precision);
			for (int j=0; j<tsa_ineq_misfit.length; j++)
				assertEquals("solution doesn't match at index "+j, sa_ineq_misfit[j], tsa_ineq_misfit[j], precision);
		}
	}

}
