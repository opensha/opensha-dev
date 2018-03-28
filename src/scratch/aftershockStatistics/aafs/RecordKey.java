package scratch.aftershockStatistics.aafs;

import org.bson.types.ObjectId;

/**
 * Key for establishing a correspondence between records in different database tables.
 * Author: Michael Barall 03/18/2018.
 *
 * In the AAFS server there are cases where there is a correspondence between
 * records in different tables (MongoDB collections), for example, an aftershock
 * forecast in one table and a cached aftershock sequence in another table.
 * If you possess one such record, you can obtain its RecordKey and then use
 * it to query the other table.
 *
 * Only code very close to the database engine should create these objects or
 * access their contents.  All other code should treat these objects as opaque.
 * (If this was C++, everything would be private and the classes that need access
 * would be friends, but Java does not have an analogous mechanism.)
 */
public class RecordKey {

	// The ObjectId, which can be used as a key for MongoDB.

    private ObjectId id;

	// Constructor saves the ObjectId.

	public RecordKey (ObjectId id) {
		this.id = id;
	}

	// Get the ObjectId.

	public ObjectId getId () {
		return id;
	}




	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "RecordKey: " + ((id == null) ? ("null") : (id.toHexString()));
		return str;
	}




	// create_unique_key - Create a key, whose value is different than any other key.

	public static RecordKey create_unique_key () {
		return new RecordKey (new ObjectId());
	}
}
