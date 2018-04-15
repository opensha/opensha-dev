package scratch.aftershockStatistics.aafs.entity;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import org.bson.types.ObjectId;
import org.mongodb.morphia.annotations.Entity;
import org.mongodb.morphia.annotations.Id;
import org.mongodb.morphia.annotations.Index;
import org.mongodb.morphia.annotations.Indexed;
import org.mongodb.morphia.annotations.IndexOptions;
import org.mongodb.morphia.annotations.Indexes;
import org.mongodb.morphia.annotations.Field;
import org.mongodb.morphia.annotations.Transient;
import org.mongodb.morphia.annotations.Embedded;

import org.mongodb.morphia.Datastore;
import org.mongodb.morphia.FindAndModifyOptions;

import org.mongodb.morphia.query.Query;
import org.mongodb.morphia.query.UpdateOperations;
import org.mongodb.morphia.query.MorphiaIterator;

import scratch.aftershockStatistics.aafs.DBCorruptException;
import scratch.aftershockStatistics.aafs.MongoDBUtil;
import scratch.aftershockStatistics.aafs.RecordKey;
import scratch.aftershockStatistics.aafs.RecordPayload;
import scratch.aftershockStatistics.aafs.RecordIterator;

import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.CompactEqkRupList;




/**
 * Holds an AAFS earthquake catalog snapshot in the MongoDB database.
 * Author: Michael Barall 04/02/2018.
 *
 * The collection "catalog" holds the earthquake catalog snapshots.
 */
@Entity(value = "catalog", noClassnameStored = true)
public class CatalogSnapshot implements java.io.Serializable {

	//----- Envelope information -----

	// Globally unique identifier for this catalog snapshot.
	// This is the MongoDB identifier.
	// Note that ObjectId implements java.io.Serializable.
	// This is set to the same value as the id of the task that generated the snapshot.

    @Id
    private ObjectId id;

	//----- Catalog information -----

	// Event ID that this earthquake sequence pertains to.

	private String event_id;

	// Start time of this earthquake sequence, in milliseconds since the epoch.
	// This is not the time of the first earthquake, because there may be an interval of no earthquakes.

	private long start_time;

	// End time of this earthquake sequence, in milliseconds since the epoch.
	// This is not the time of the last earthquake, because there may be an interval of no earthquakes.

	private long end_time;

	// eqk_count - Number of earthquakes.

	private int eqk_count;

	// lat_lon_depth_list - Compressed latitude, longitude, and depth for each earthquake.
	// An entry may be zero if there is no location data for the particular earthquake.
	// The length must equal max(eqk_count,1).

	@Embedded
	private long[] lat_lon_depth_list;

	// mag_time_list - Compressed magnitude and time for each earthquake.
	// An entry may be zero if there is no time and magnitude data for the particular earthquake.
	// The length must equal max(eqk_count,1).

	@Embedded
	private long[] mag_time_list;




	//----- Getters and setters -----

    private ObjectId get_id() {
        return id;
    }

    private void set_id (ObjectId id) {
        this.id = id;
    }

    public String get_event_id() {
        return event_id;
    }

    private void set_event_id (String event_id) {
        this.event_id = event_id;
    }

    public long get_start_time() {
        return start_time;
    }

    private void set_start_time (long start_time) {
        this.start_time = start_time;
    }

    public long get_end_time() {
        return end_time;
    }

    private void set_end_time (long end_time) {
        this.end_time = end_time;
    }

    public int get_eqk_count() {
        return eqk_count;
    }

    private void set_eqk_count (int eqk_count) {
        this.eqk_count = eqk_count;
    }

    private long[] get_lat_lon_depth_list() {
        return lat_lon_depth_list;
    }

    private void set_lat_lon_depth_list (long[] lat_lon_depth_list) {
        this.lat_lon_depth_list = lat_lon_depth_list;
    }

    private long[] get_mag_time_list() {
        return mag_time_list;
    }

    private void set_mag_time_list (long[] mag_time_list) {
        this.mag_time_list = mag_time_list;
    }




	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "CatalogSnapshot\n"
			+ "\tid: " + ((id == null) ? ("null") : (id.toHexString())) + "\n"
			+ "\tevent_id: " + event_id + "\n"
			+ "\tstart_time: " + start_time + "\n"
			+ "\tend_time: " + end_time + "\n"
			+ "\teqk_count: " + eqk_count + "\n"
			+ "\tlat_lon_depth_list: " + ((lat_lon_depth_list == null) ? ("null") : ("len=" + lat_lon_depth_list.length)) + "\n"
			+ "\tmag_time_list: " + ((mag_time_list == null) ? ("null") : ("len=" + mag_time_list.length));
		return str;
	}




	/**
	 * get_record_key - Get the record key for this catalog snapshot.
	 */
	public RecordKey get_record_key () {
		return new RecordKey(id);
	}




	/**
	 * set_record_key - Set the record key for this catalog snapshot.
	 * @param key = Record key. Can be null.
	 */
	private void set_record_key (RecordKey key) {
		if (key == null) {
			id = null;
		} else {
			id = key.getId();
		}
		return;
	}




	/**
	 * get_rupture_list - Get the earthquake rupture list for this catalog snapshot.
	 */
	public CompactEqkRupList get_rupture_list () {

		// For empty list, pass in zero-size arrays

		if (eqk_count == 0) {
			return new CompactEqkRupList (eqk_count, new long[0], new long[0]);
		}

		// For non-empty list, pass our arrays

		return new CompactEqkRupList (eqk_count, lat_lon_depth_list, mag_time_list);
	}




	/**
	 * set_rupture_list - Set the earthquake rupture list for this catalog snapshot.
	 * @param rupture_list = Rupture list. Cannot be null.
	 */
	private void set_rupture_list (CompactEqkRupList rupture_list) {
		eqk_count = rupture_list.get_eqk_count();

		// For empty list, use one-element lists

		if (eqk_count == 0) {
			lat_lon_depth_list = new long[1];
			lat_lon_depth_list[0] = 0L;
			mag_time_list = new long[1];
			mag_time_list[0] = 0L;
			return;
		}

		// For non-empty list, pull the existing arrays, and re-size them if needed

		lat_lon_depth_list = rupture_list.get_lat_lon_depth_list();
		mag_time_list = rupture_list.get_mag_time_list();

		if (lat_lon_depth_list.length != eqk_count) {
			lat_lon_depth_list = Arrays.copyOf (lat_lon_depth_list, eqk_count);
		}
		if (mag_time_list.length != eqk_count) {
			mag_time_list = Arrays.copyOf (mag_time_list, eqk_count);
		}
		return;
	}




	/**
	 * submit_catalog_shapshot - Submit a catalog snapshot.
	 * @param key = Record key associated with this catalog snapshot. Can be null to assign a new one.
	 * @param event_id = Event associated with this catalog snapshot, or "" if none. Cannot be null.
	 * @param start_time = Start time of this earthquake sequence, in milliseconds
	 *                     since the epoch. Must be positive.
	 * @param end_time = End time of this earthquake sequence, in milliseconds
	 *                   since the epoch. Must be positive. Must be >= start_time.
	 * @param rupture_list = Rupture list. Cannot be null.
	 * @return
	 * Returns the new entry.
	 */
	public static CatalogSnapshot submit_catalog_shapshot (RecordKey key, String event_id,
			long start_time, long end_time, CompactEqkRupList rupture_list) {

		// Check conditions

		if (!( event_id != null
			&& start_time > 0L
			&& end_time >= start_time
			&& rupture_list != null )) {
			throw new IllegalArgumentException("CatalogSnapshot.submit_catalog_shapshot: Invalid catalog snapshot parameters");
		}

		// Construct the catalog snapshot object

		CatalogSnapshot catsnap = new CatalogSnapshot();
		catsnap.set_record_key (key);
		catsnap.set_event_id (event_id);
		catsnap.set_start_time (start_time);
		catsnap.set_end_time (end_time);
		catsnap.set_rupture_list (rupture_list);

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(catsnap);
		
		return catsnap;
	}




	/**
	 * store_catalog_shapshot - Store a catalog snapshot into the database.
	 * This is primarily for restoring from backup.
	 */
	public static CatalogSnapshot store_catalog_shapshot (CatalogSnapshot catsnap) {

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(catsnap);
		
		return catsnap;
	}




	/**
	 * get_catalog_shapshot_for_key - Get the catalog snapshot with the given key.
	 * @param key = Record key. Cannot be null or empty.
	 * Returns the catalog snapshot, or null if not found.
	 */
	public static CatalogSnapshot get_catalog_shapshot_for_key (RecordKey key) {

		if (!( key != null && key.getId() != null )) {
			throw new IllegalArgumentException("CatalogSnapshot.get_catalog_shapshot_for_key: Missing or empty record key");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select id == key

		Query<CatalogSnapshot> query = datastore.createQuery(CatalogSnapshot.class)
											.filter("id ==", key.getId());

		// Run the query

		CatalogSnapshot catsnap = query.get();

		return catsnap;
	}




	/**
	 * get_catalog_snapshot_range - Get a range of catalog snapshots, reverse-sorted by end time.
	 * @param end_time_lo = Minimum end time, in milliseconds since the epoch.
	 *                      Can be 0L for no minimum.
	 * @param end_time_hi = Maximum end time, in milliseconds since the epoch.
	 *                      Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 */
	public static List<CatalogSnapshot> get_catalog_snapshot_range (long end_time_lo, long end_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<CatalogSnapshot> query = datastore.createQuery(CatalogSnapshot.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with end_time >= end_time_lo

		if (end_time_lo > 0L) {
			query = query.filter("end_time >=", new Long(end_time_lo));
		}

		// Select entries with end_time <= end_time_hi

		if (end_time_hi > 0L) {
			query = query.filter("end_time <=", new Long(end_time_hi));
		}

		// Sort by end_time in descending order (most recent first)

		query = query.order("-end_time");

		// Run the query

		List<CatalogSnapshot> entries = query.asList();

		return entries;
	}




	/**
	 * fetch_catalog_snapshot_range - Iterate a range of catalog snapshots, reverse-sorted by end time.
	 * @param end_time_lo = Minimum end time, in milliseconds since the epoch.
	 *                      Can be 0L for no minimum.
	 * @param end_time_hi = Maximum end time, in milliseconds since the epoch.
	 *                      Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 */
	public static RecordIterator<CatalogSnapshot> fetch_catalog_snapshot_range (long end_time_lo, long end_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<CatalogSnapshot> query = datastore.createQuery(CatalogSnapshot.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with end_time >= end_time_lo

		if (end_time_lo > 0L) {
			query = query.filter("end_time >=", new Long(end_time_lo));
		}

		// Select entries with end_time <= end_time_hi

		if (end_time_hi > 0L) {
			query = query.filter("end_time <=", new Long(end_time_hi));
		}

		// Sort by end_time in descending order (most recent first)

		query = query.order("-end_time");

		// Run the query

		MorphiaIterator<CatalogSnapshot, CatalogSnapshot> morphia_iterator = query.fetch();

		return new RecordIterator<CatalogSnapshot>(morphia_iterator);
	}




	/**
	 * delete_catalog_snapshot - Delete a catalog snapshot.
	 * @param entry = Existing catalog snapshot to delete.
	 * @return
	 */
	public static void delete_catalog_snapshot (CatalogSnapshot entry) {

		// Check conditions

		if (!( entry != null && entry.get_id() != null )) {
			throw new IllegalArgumentException("CatalogSnapshot.delete_catalog_snapshot: Invalid parameters");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Run the delete

		datastore.delete(entry);
		
		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 10001;

	private static final String M_VERSION_NAME = "CatalogSnapshot";

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		String sid = id.toHexString();
		writer.marshalString      ("id"                , sid               );
		writer.marshalString      ("event_id"          , event_id          );
		writer.marshalLong        ("start_time"        , start_time        );
		writer.marshalLong        ("end_time"          , end_time          );
		writer.marshalInt         ("eqk_count"         , eqk_count         );
		writer.marshalLongArray   ("lat_lon_depth_list", lat_lon_depth_list);
		writer.marshalLongArray   ("mag_time_list"     , mag_time_list     );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		String sid;
		sid                = reader.unmarshalString      ("id"                );
		event_id           = reader.unmarshalString      ("event_id"          );
		start_time         = reader.unmarshalLong        ("start_time"        );
		end_time           = reader.unmarshalLong        ("end_time"          );
		eqk_count          = reader.unmarshalInt         ("eqk_count"         );
		lat_lon_depth_list = reader.unmarshalLongArray   ("lat_lon_depth_list");
		mag_time_list      = reader.unmarshalLongArray   ("mag_time_list"     );
		id = new ObjectId(sid);

		return;
	}

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public CatalogSnapshot unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}


}
