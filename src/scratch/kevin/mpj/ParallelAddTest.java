package scratch.kevin.mpj;

import mpi.Datatype;
import mpi.MPI;
import mpi.MPIException;

public class ParallelAddTest {
	
	public static void main(String[] args) {
		try {
			MPI.Init(args);
			
			int me = MPI.COMM_WORLD.Rank();
			int size = MPI.COMM_WORLD.Size();
			
//		MPI.COMM_WORLD.Bsend(buf, offset, count, datatype, dest, tag)
			int offset = 0;
			int count = 200;
			Datatype datatype = MPI.DOUBLE;
			int dest = 0;
			int tag = 0;
			if (me > 0) {
				double[] buf = new double[count];
				for (int i=0; i<count; i++) {
					buf[i] = me*1000 + count;
				}
				MPI.COMM_WORLD.Send(buf, offset, count, datatype, dest, tag);
			} else {
				int cnt = 0;
				int should = 0;
				for (int proc=1; proc<size; proc++) {
					double[] buf = new double[count];
					MPI.COMM_WORLD.Recv(buf, offset, count, datatype, proc, tag);
					for (double val : buf)
						cnt += val;
					should += proc;
				}
				System.out.println("Total count: "+cnt);
//				System.out.println("Should be: "+should);
			}
			
			MPI.Finalize();
		} catch (MPIException e) {
			MPI.COMM_WORLD.Abort(1);
		}
	}

}
