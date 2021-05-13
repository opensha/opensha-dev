package scratch.kevin.ucerf3;

import java.sql.SQLException;

import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionConnectionsDB_DAO;
import org.opensha.refFaultParamDb.dao.exception.QueryException;
import org.opensha.refFaultParamDb.vo.FaultSectionConnection;
import org.opensha.refFaultParamDb.vo.FaultSectionConnectionList;

public class ConnectionDBTest {
	
	private static void printConnections(FaultSectionConnectionList conns) {
		for (FaultSectionConnection conn : conns) {
			System.out.println(conn.toString());
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		DB_AccessAPI db = null;
		
		try {
			db= DB_ConnectionPool.getDirectLatestReadWriteConnection();
			
			FaultSectionConnectionsDB_DAO connsDB = new FaultSectionConnectionsDB_DAO(db);
			
			System.out.println("Initial connections");
			printConnections(connsDB.getAllConnections());
			
			System.out.println("Adding connection.");
			FaultSectionConnection conn = new FaultSectionConnection(0, 1, new Location(34, -118), new Location(34.1, -118));
			connsDB.addConnection(conn);
			
			System.out.println("Printing connections");
			printConnections(connsDB.getAllConnections());
			
//			System.out.println("Removing connections");
////			connsDB.removeAllConnections(0);
//			connsDB.removeConnection(1, 0);
//			
//			System.out.println("Printing connections");
//			printConnections(connsDB.getAllConnections());
			
			System.out.println("DONE!");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if (db != null) {
				try {
					db.destroy();
					db = null;
				} catch (SQLException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		System.exit(0);
	}

}
