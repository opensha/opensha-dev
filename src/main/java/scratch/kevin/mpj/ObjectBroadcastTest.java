package scratch.kevin.mpj;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import mpi.MPI;
import mpi.MPIException;

public class ObjectBroadcastTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			MPI.Init(args);
			
			int me = MPI.COMM_WORLD.Rank();
			int size = MPI.COMM_WORLD.Size();
			
			Site[] buf = new Site[1];
			
			if (me == 0) {
				Site site = new Site(new Location(34, -118));
				Vs30_Param vs30 = new Vs30_Param(760);
				vs30.setValueAsDefault();
				site.addParameter(vs30);
				buf[0] = site;
			} else {
				buf[0] = null;
			}
			
			MPI.COMM_WORLD.Bcast(buf, 0, 1, MPI.OBJECT, 0);
			
			System.out.println("Rank: "+me+"\n"+buf[0]);
			
			MPI.Finalize();
		} catch (MPIException e) {
			e.printStackTrace();
			MPI.COMM_WORLD.Abort(1);
		}
	}

}
