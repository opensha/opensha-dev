package scratch.kevin.nshm23.dmWeightInversion;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.opensha.commons.calc.nnls.NNLSWrapper;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.Quantity;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfits;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class TestInversion {

	public static void main(String[] args) throws IOException {
		List<double[]> dms = new ArrayList<>();
		List<String> names = new ArrayList<>();
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v2;
		double[] median = null;
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			if (dm == NSHM23_DeformationModels.MEDIAN) {
				List<? extends FaultSection> sects = dm.build(fm);
				median = new double[sects.size()];
				for (int s=0; s<sects.size(); s++)
					median[s] = sects.get(s).getOrigAveSlipRate();
			} else if (dm.getNodeWeight(null) > 0) {
				List<? extends FaultSection> sects = dm.build(fm);
				double[] slips = new double[sects.size()];
				for (int s=0; s<sects.size(); s++)
					slips[s] = sects.get(s).getOrigAveSlipRate();
				dms.add(slips);
				names.add(dm.getName());
			}
		}
		double[] initial = new double[dms.size()];
		for (int i=0; i<initial.length; i++)
			initial[i] = 1d/initial.length;
		
		List<InversionConstraint> constraints = new ArrayList<>();
		constraints.add(new MatchMedianConstraint(dms, median, 0.1, ConstraintWeightingType.NORMALIZED));
		constraints.add(new MatchMedianConstraint(dms, median, 1, ConstraintWeightingType.UNNORMALIZED));
		constraints.add(new SumToOneConstraint(10000));
		
		int rows = 0;
		for (InversionConstraint constraint : constraints)
			rows += constraint.getNumRows();
		DoubleMatrix2D A = new SparseDoubleMatrix2D(rows, dms.size());
		double[] d = new double[A.rows()];
		
		int row = 0;
		List<ConstraintRange> ranges = new ArrayList<>();
		for (InversionConstraint constraint : constraints) {
			constraint.encode(A, d, row);
			ranges.add(constraint.getRange(row));
			row += constraint.getNumRows();
		}
		
		SerialSimulatedAnnealing sa = new SerialSimulatedAnnealing(A, d, initial);
		
		sa.setPerturbationFunc(GenerationFunctionType.UNIFORM_0p001);
		System.out.println("Iterating...");
		sa.iterate(1000000l);
		
		double[] sol = sa.getBestSolution();
		
		InversionMisfits misfits = new InversionMisfits(ranges, sa.getBestMisfit(), d);
		
		for (MisfitStats stats : misfits.getMisfitStats().getStats()) {
			System.out.println(stats.range.name+" Misfits:");
			for (Quantity q : Quantity.values())
				System.out.println("\t"+q+":\t"+(float)stats.get(q));
		}
		
		System.out.println("Final weights:");
		double sum = 0d;
		for (int i=0; i<sol.length; i++) {
			System.out.println(names.get(i)+": "+(float)sol[i]);
			sum += sol[i];
		}
		System.out.println("SUM: "+sum);
		
		NNLSWrapper nnls = new NNLSWrapper();

		int nRow = A.rows();
		int nCol = A.columns();
		
//		System.out.println("NNLS: nRow="+nRow+"; nCol="+nCol);
		
		double[] Amat = new double[nRow*nCol];
		double[] x = new double[nCol];
			
		int k = 0;
		for(int j=0;j<nCol;j++) 
			for(int i=0; i<nRow;i++)	{
				Amat[k]=A.get(i, j);
				k+=1;
			}
		nnls.update(Amat,nRow,nCol);
		
		boolean converged = nnls.solve(d,x);
		if(!converged)
			throw new RuntimeException("ERROR:  NNLS Inversion Failed");
		
//		System.out.println("A value: "+x[x.length-1]);
		System.out.println("NNLS weights:");
		sum = 0d;
		for (int i=0; i<sol.length; i++) {
			System.out.println(names.get(i)+": "+(float)x[i]);
			sum += x[i];
		}
		System.out.println("SUM: "+sum);
	}

}
