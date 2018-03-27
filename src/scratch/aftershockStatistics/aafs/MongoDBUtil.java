package scratch.aftershockStatistics.aafs;

import org.mongodb.morphia.Datastore;
import org.mongodb.morphia.Morphia;
import org.mongodb.morphia.mapping.MapperOptions;
import org.mongodb.morphia.mapping.Mapper;

import com.mongodb.MongoClient;
import com.mongodb.ServerAddress;
import com.mongodb.MongoCredential;
import com.mongodb.MongoClientOptions;
import com.mongodb.client.MongoDatabase;


/**
 * This class holds utilities for access to MongoDB.
 * Author: Michael Barall 03/15/2018.
 *
 * Any operation that uses MongoDB must create one of these objects.
 * It is strongly advised to create the object in a try-with-resources
 * statement, to ensure it is closed upon exit from the program.
 */
public class MongoDBUtil implements AutoCloseable {

	// The MongoDB server address.

	private static ServerAddress serverAddress = null;

	// The MongoDB security credentials, used to log in to the MongoDB server.

	private static MongoCredential credentials = null;

	// The MongoDB client options.

	private static MongoClientOptions mongoOptions = null;

	// The MongoDB client endpoint.
	// Note that there should be only one client per JVM.

	private static MongoClient mongoClient = null;

	// The MongoDB database.

    //private static MongoDatabase db = null;

	// The Morphia endpoint.

	private static Morphia morphia = null;

	// The Morphia datastore.

    private static Datastore datastore = null;



	
	/**
	 * Attach to the MongoDB database.
	 */
	public MongoDBUtil() {

		// This cannot be called if the connection is currently open.

		if (mongoClient != null) {
			throw new RuntimeException("MongoDBUtil: Connection to MongoDB is already open");
		}

		MongoClient saved_mongoClient = null;

		try {

			// Get the server configuration, which has the database address and credentials.

			ServerConfig config = new ServerConfig();

			// Create the address of the server, using host IP address and port.
			// Note: ServerAddress offers several ways to specify the address.

			serverAddress = new ServerAddress(config.getDb_host(), config.getDb_port());

			// Create the login credentials, for username, database name, and password.
			// Note: MongoCredential can create various other sorts of credentials.
			// Note: In MongoDB, it is necessary to authenticate to a particular database.
			//  It must be the database that was used to create the user account.
			//  This does not limit the databases that can be used, once logged in.

			credentials = MongoCredential.createCredential(config.getDb_user(), config.getDb_name(), config.getDb_password().toCharArray());

			// Create the MongoDB client options.
			// Note: We use the default client options.
			// Note: MongoClientOptions offers many options that can be set, using manipulator methods.
			// For example, this would set the connection timeout to connectTimeout milliseconds:
			//  new MongoClientOptions.Builder().connectTimeout(connectTimeout).build();

			mongoOptions = new MongoClientOptions.Builder().build();

			// Create the MongoDB client endpoint.

			mongoClient = new MongoClient(serverAddress, credentials, mongoOptions);
			saved_mongoClient = mongoClient;

			// Apparently the Mongo client lazy connects, this call forces check for connection success.

			mongoClient.getAddress();

			// Get the database, using database name.
			// This could be used for database operations not supported by Morphia.

			//db = mongoClient.getDatabase(config.getDb_name());

			// Create the Morphia endpoint.

			morphia = new Morphia();

			// At this point we could configure mapping options.
			// The most common options to be configured are storeEmpties (which selects whether
			// empty List, Map, Set, and array values are stored; default false), and storeNulls
			// (which selects whether null values are stored; default false).
			// Apparently this would be done like so:
			//  MapperOptions options = new MapperOptions();
			//  options.setStoreEmpties(true);
			//  options.setStoreNulls(true);
			//  morphia.getMapper().setOptions(options);

			// Tell Morphia where to find our classes.
			// Morphia finds every class in the specified package that is annotated with @Entity,
			// and reads its metadata.
			// This could be called multiple times to specify multiple packages.

			morphia.mapPackage("scratch.aftershockStatistics.aafs.entity");

			// Create the Morphia datastore, using the database name.

			datastore = morphia.createDatastore(mongoClient, config.getDb_name());

			// This ensures the existence of any indexes found during class mapping.
			// Indexes are created if necessary.

			datastore.ensureIndexes();

		} catch (Exception e) {
			datastore = null;
			morphia = null;
			//db = null;
			mongoClient = null;
			if (saved_mongoClient != null) {
				try {
					saved_mongoClient.close();
				} catch (Exception e2) {
				}
			}
			saved_mongoClient = null;
			mongoOptions = null;
			credentials = null;
			serverAddress = null;
			throw new RuntimeException("MongoDBUtil: Unable to connect to MongoDB", e);
		}
	}



	
	/**
	 * Close the MongoDB database.
	 */
	@Override
	public void close() {

		// Close the client if it is currently open.

		if (mongoClient != null) {
			datastore = null;
			morphia = null;
			//db = null;
			mongoClient.close();
			mongoClient = null;
			mongoOptions = null;
			credentials = null;
			serverAddress = null;
		}

		return;
	}




	
	/**
	 * Retrieve the MongoDB database.
	 */
    //public static DB getDB() {
    //    return db;
    //}

 


	
	/**
	 * Retrieve the Morphia datastore.
	 */
   public static Datastore getDatastore() {
        return datastore;
    }
}
