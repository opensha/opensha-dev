package scratch.aftershockStatistics.aafs.entity;

import java.util.List;

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




/**
 * Holds an AAFS timeline entry in the MongoDB database.
 * Author: Michael Barall 04/08/2018.
 *
 * The collection "timeline" holds the timeline entries.
 */
@Entity(value = "timeline", noClassnameStored = true)
@Indexes({
	@Index(fields = {@Field("event_id")}, options = @IndexOptions(name = "actevid")),
	@Index(fields = {@Field("action_time")}, options = @IndexOptions(name = "acttime"))
})
public class TimelineEntry implements java.io.Serializable {

	//----- Envelope information -----

	// Globally unique identifier for this timeline entry.
	// This is the MongoDB identifier.
	// Note that ObjectId implements java.io.Serializable.
	// This is set to the same value as the id of the task that generated the timeline entry.

    @Id
    private ObjectId id;

	// Time that this timeline entry was created, in milliseconds since the epoch.
	// The collection is indexed on this field, so timeline entries in a given
	// time span can be obtained with a database query efficiently.

	//@Indexed(options = @IndexOptions(name = "acttime"))
	private long action_time;

	//----- Action information -----

	// Event ID that this action pertains to.
	// Entries not referring to an event should put an empty string here (not null).
	// This can be used to find all timeline entries pertaining to an event.

	//@Indexed(options = @IndexOptions(name = "actevid"))
	private String event_id;

	// Action code for this timeline entry.

	private int actcode;

	// Details of this action.
	// Any additional information needed is stored as a JSON string containing marshaled data.
	// If none, this should be an empty string (not null).

	private String details;

//	// Details of this action.
//	// Any additional information needed is stored as marshaled data.
//	// Each array should have at least one element.
//
//	@Embedded
//	private long[] details_l;
//	@Embedded
//	private double[] details_d;
//	@Embedded
//	private String[] details_s;




	//----- Getters and setters -----

    private ObjectId get_id() {
        return id;
    }

    private void set_id (ObjectId id) {
        this.id = id;
    }

    public long get_action_time() {
        return action_time;
    }

    private void set_action_time (long action_time) {
        this.action_time = action_time;
    }

    public String get_event_id() {
        return event_id;
    }

    private void set_event_id (String event_id) {
        this.event_id = event_id;
    }

    public int get_actcode() {
        return actcode;
    }

    private void set_actcode (int actcode) {
        this.actcode = actcode;
    }




	/**
	 * get_details - Get a reader for the details.
	 */
    public MarshalReader get_details() {
		Object json_source;
		if (details == null) {
			json_source = null;
		}
		else if (details.equals("")) {
			json_source = null;
		}
		else {
			json_source = details;
		}
        return new MarshalImpJsonReader (json_source);
    }


	/**
	 * set_details - Set details from the marshaled data, can be null for none.
	 */
    private void set_details (MarshalWriter writer) {

		if (writer == null) {
			details = "";
			return;
		}

		if (writer instanceof MarshalImpJsonWriter) {
			MarshalImpJsonWriter w = (MarshalImpJsonWriter)writer;
			if (w.check_write_complete()) {
				details = w.get_json_string();
			} else {
				details = "";
			}
			return;
		}

		if (writer instanceof RecordPayload) {
			RecordPayload p = (RecordPayload)writer;
			details = p.get_json_string();
			return;
		}

		throw new IllegalArgumentException("TimelineEntry.set_details: Incorrect type of marshal writer");
    }


	/**
	 * begin_details - Get a writer to use for marhaling details.
	 */
    public static MarshalWriter begin_details() {
        return new MarshalImpJsonWriter ();
    }


	/**
	 * get_details_as_payload - Get a writer containing the details.
	 */
    RecordPayload get_details_as_payload() {
        return new RecordPayload (details);
    }


	/**
	 * get_details_description - Get a string describing the details.
	 */
	private String get_details_description () {
		return ((details == null) ? "null" : ("len = " + details.length()));
	}


	/**
	 * dump_details - Dump details into a string, for trouble-shooting.
	 */
	public String dump_details () {
		return ((details == null) ? "null" : details);
	}




//	/**
//	 * get_details - Get a reader for the details.
//	 */
//    public MarshalReader get_details() {
//        return new MarshalImpArray (details_l, details_d, details_s);
//    }
//
//
//	/**
//	 * set_details - Set details from the marshaled data, can be null for none.
//	 */
//    private void set_details (MarshalWriter writer) {
//
//		if (writer == null) {
//			details_l = new long[1];
//			details_l[0] = 0L;
//			details_d = new double[1];
//			details_d[0] = 0.0;
//			details_s = new String[1];
//			details_s[0] = "";
//			return;
//		}
//
//		if (writer instanceof MarshalImpArray) {
//			MarshalImpArray w = (MarshalImpArray)writer;
//			if (w.check_write_complete()) {
//				details_l = w.get_long_store();
//				details_d = w.get_double_store();
//				details_s = w.get_string_store();
//			} else {
//				details_l = new long[1];
//				details_l[0] = 0L;
//				details_d = new double[1];
//				details_d[0] = 0.0;
//				details_s = new String[1];
//				details_s[0] = "";
//			}
//			return;
//		}
//
//		if (writer instanceof RecordPayload) {
//			RecordPayload p = (RecordPayload)writer;
//			details_l = p.get_long_store();
//			details_d = p.get_double_store();
//			details_s = p.get_string_store();
//			return;
//		}
//
//		throw new IllegalArgumentException("PendingTask.set_details: Incorrect type of marshal writer");
//    }
//
//
//	/**
//	 * begin_details - Get a writer to use for marhaling details.
//	 */
//    public static MarshalWriter begin_details() {
//        return new MarshalImpArray ();
//    }
//
//
//	/**
//	 * get_details_as_payload - Get a writer containing the details.
//	 */
//    RecordPayload get_details_as_payload() {
//        return new RecordPayload (details_l, details_d, details_s);
//    }
//
//
//	/**
//	 * get_details_description - Get a string describing the details.
//	 */
//	private String get_details_description () {
//		return "llen = " + details_l.length + ", dlen = " + details_d.length + ", slen = " + details_s.length;
//	}
//
//
//	///**
//	// * These getters are for use by other package members for direct copy of the details.
//	// */
//    //long[] get_details_l() {
//    //    return details_l;
//    //}
//    //double[] get_details_d() {
//    //    return details_d;
//    //}
//    //String[] get_details_s() {
//    //    return details_s;
//    //}




	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "TimelineEntry\n"
			+ "\tid: " + ((id == null) ? ("null") : (id.toHexString())) + "\n"
			+ "\taction_time: " + action_time + "\n"
			+ "\tevent_id: " + event_id + "\n"
			+ "\tactcode: " + actcode + "\n"
			+ "\tdetails: " + get_details_description();
		return str;
	}




	/**
	 * get_record_key - Get the record key for this timeline entry.
	 */
	public RecordKey get_record_key () {
		return new RecordKey(id);
	}




	/**
	 * set_record_key - Set the record key for this timeline entry.
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
	 * submit_timeline_entry - Submit a timeline entry.
	 * @param key = Record key associated with this task. Can be null to assign a new one.
	 * @param action_time = Time of this timeline entry, in milliseconds
	 *                   since the epoch. Must be positive.
	 * @param event_id = Event associated with this task, or "" if none. Cannot be null.
	 * @param actcode = Operation code used to dispatch the task.
	 * @param details = Further details of this task. Can be null if there are none.
	 * @return
	 * Returns the new entry.
	 */
	public static TimelineEntry submit_timeline_entry (RecordKey key, long action_time, String event_id,
			int actcode, MarshalWriter details) {

		// Check conditions

		if (!( action_time > 0L
			&& event_id != null )) {
			throw new IllegalArgumentException("TimelineEntry.submit_timeline_entry: Invalid log parameters");
		}

		// Construct the timeline entry object

		TimelineEntry tentry = new TimelineEntry();
		tentry.set_record_key (key);
		tentry.set_action_time (action_time);
		tentry.set_event_id (event_id);
		tentry.set_actcode (actcode);
		tentry.set_details (details);

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(tentry);
		
		return tentry;
	}




	/**
	 * store_timeline_entry - Store a timeline entry into the database.
	 * This is primarily for restoring from backup.
	 */
	public static TimelineEntry store_timeline_entry (TimelineEntry tentry) {

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(tentry);
		
		return tentry;
	}




	/**
	 * get_timeline_entry_for_key - Get the timeline entry with the given key.
	 * @param key = Record key. Cannot be null or empty.
	 * Returns the timeline entry, or null if not found.
	 */
	public static TimelineEntry get_timeline_entry_for_key (RecordKey key) {

		if (!( key != null && key.getId() != null )) {
			throw new IllegalArgumentException("TimelineEntry.get_timeline_entry_for_key: Missing or empty record key");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select id == key

		Query<TimelineEntry> query = datastore.createQuery(TimelineEntry.class)
											.filter("id ==", key.getId());

		// Run the query

		TimelineEntry tentry = query.get();

		return tentry;
	}




	/**
	 * get_timeline_entry_range - Get a range of timeline entries, reverse-sorted by action time.
	 * @param action_time_lo = Minimum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no minimum.
	 * @param action_time_hi = Maximum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 */
	public static List<TimelineEntry> get_timeline_entry_range (long action_time_lo, long action_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<TimelineEntry> query = datastore.createQuery(TimelineEntry.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with action_time >= action_time_lo

		if (action_time_lo > 0L) {
			query = query.filter("action_time >=", new Long(action_time_lo));
		}

		// Select entries with action_time <= action_time_hi

		if (action_time_hi > 0L) {
			query = query.filter("action_time <=", new Long(action_time_hi));
		}

		// Sort by action_time in descending order (most recent first)

		query = query.order("-action_time");

		// Run the query

		List<TimelineEntry> entries = query.asList();

		return entries;
	}




	/**
	 * fetch_timeline_entry_range - Iterate a range of timeline entries, reverse-sorted by action time.
	 * @param action_time_lo = Minimum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no minimum.
	 * @param action_time_hi = Maximum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 */
	public static RecordIterator<TimelineEntry> fetch_timeline_entry_range (long action_time_lo, long action_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<TimelineEntry> query = datastore.createQuery(TimelineEntry.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with action_time >= action_time_lo

		if (action_time_lo > 0L) {
			query = query.filter("action_time >=", new Long(action_time_lo));
		}

		// Select entries with action_time <= action_time_hi

		if (action_time_hi > 0L) {
			query = query.filter("action_time <=", new Long(action_time_hi));
		}

		// Sort by action_time in descending order (most recent first)

		query = query.order("-action_time");

		// Run the query

		MorphiaIterator<TimelineEntry, TimelineEntry> morphia_iterator = query.fetch();

		return new RecordIterator<TimelineEntry>(morphia_iterator);
	}




	/**
	 * get_recent_timeline_entry - Get the most recent in a range of timeline entries.
	 * @param action_time_lo = Minimum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no minimum.
	 * @param action_time_hi = Maximum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 * Returns the matching timeline entry with the greatest action_time (most recent),
	 * or null if there is no matching timeline entry.
	 */
	public static TimelineEntry get_recent_timeline_entry (long action_time_lo, long action_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<TimelineEntry> query = datastore.createQuery(TimelineEntry.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with action_time >= action_time_lo

		if (action_time_lo > 0L) {
			query = query.filter("action_time >=", new Long(action_time_lo));
		}

		// Select entries with action_time <= action_time_hi

		if (action_time_hi > 0L) {
			query = query.filter("action_time <=", new Long(action_time_hi));
		}

		// Sort by action_time in descending order (most recent first)

		query = query.order("-action_time");

		// Run the query

		TimelineEntry tentry = query.get();

		return tentry;
	}




	/**
	 * delete_timeline_entry - Delete a timeline entry.
	 * @param entry = Existing timeline entry to delete.
	 * @return
	 */
	public static void delete_timeline_entry (TimelineEntry entry) {

		// Check conditions

		if (!( entry != null && entry.get_id() != null )) {
			throw new IllegalArgumentException("TimelineEntry.delete_timeline_entry: Invalid parameters");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Run the delete

		datastore.delete(entry);
		
		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 11001;

	private static final String M_VERSION_NAME = "TimelineEntry";

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		String sid = id.toHexString();
		writer.marshalString      ("id"         , sid        );
		writer.marshalLong        ("action_time", action_time);
		writer.marshalString      ("event_id"   , event_id   );
		writer.marshalInt         ("actcode"    , actcode    );
		writer.marshalJsonString  ("details"    , details    );
//		writer.marshalLongArray   ("details_l"  , details_l  );
//		writer.marshalDoubleArray ("details_d"  , details_d  );
//		writer.marshalStringArray ("details_s"  , details_s  );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		String sid;
		sid         = reader.unmarshalString      ("id"         );
		action_time = reader.unmarshalLong        ("action_time");
		event_id    = reader.unmarshalString      ("event_id"   );
		actcode     = reader.unmarshalInt         ("actcode"    );
		details     = reader.unmarshalJsonString  ("details"    );
//		details_l   = reader.unmarshalLongArray   ("details_l"  );
//		details_d   = reader.unmarshalDoubleArray ("details_d"  );
//		details_s   = reader.unmarshalStringArray ("details_s"  );
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

	public TimelineEntry unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}


}
