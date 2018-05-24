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
 * Holds an AAFS log entry in the MongoDB database.
 * Author: Michael Barall 03/15/2018.
 *
 * The collection "log" holds the log entries.
 */
@Entity(value = "log", noClassnameStored = true)
@Indexes({
	@Index(fields = {@Field("event_id")}, options = @IndexOptions(name = "logevid")),
	@Index(fields = {@Field("log_time")}, options = @IndexOptions(name = "logtime"))
})
public class LogEntry implements java.io.Serializable {

	//----- Envelope information -----

	// Globally unique identifier for this log entry.
	// This is the MongoDB identifier.
	// Note that ObjectId implements java.io.Serializable.
	// This is set to the same value as the id of the task that generated the log entry.

    @Id
    private ObjectId id;

	// Time that this log entry was created, in milliseconds since the epoch.
	// The collection is indexed on this field, so log entries in a given
	// time span can be obtained with a database query efficiently.

	//@Indexed(options = @IndexOptions(name = "logtime"))
	private long log_time;

	//----- Task information -----

	// Event ID that this task pertains to.
	// Entries not referring to an event should put an empty string here (not null).
	// This can be used to find all log entries pertaining to an event.

	//@Indexed(options = @IndexOptions(name = "logevid"))
	private String event_id;

	// Time this task was originally scheduled to execute, in milliseconds since the epoch.

	private long sched_time;

	// Time this task was submitted, in milliseconds since the epoch.

	private long submit_time;

	// Person or entity that submitted the task.

	private String submit_id;

	// Operation code for this task.

	private int opcode;

	// Stage number, assigned by the user.
	// This can be used for a sequence of commands with the same PendingTask.

	private int stage;

	// Details of this task.
	// Any additional information needed is stored as a JSON string containing marshaled data.
	// If none, this should be an empty string (not null).

	private String details;

//	// Details of this task.
//	// Any additional information needed is stored as marshaled data.
//	// Each array should have at least one element.
//
//	@Embedded
//	private long[] details_l;
//	@Embedded
//	private double[] details_d;
//	@Embedded
//	private String[] details_s;

	//----- Result information -----

	// Result code for this task.

	private int rescode;

	// Results for this task.
	// If none, this should be an empty string (not null).

	private String results;




	//----- Getters and setters -----

    private ObjectId get_id() {
        return id;
    }

    private void set_id (ObjectId id) {
        this.id = id;
    }

    public long get_log_time() {
        return log_time;
    }

    private void set_log_time (long log_time) {
        this.log_time = log_time;
    }

    public String get_event_id() {
        return event_id;
    }

    private void set_event_id (String event_id) {
        this.event_id = event_id;
    }

    public long get_sched_time() {
        return sched_time;
    }

    private void set_sched_time (long sched_time) {
        this.sched_time = sched_time;
    }

    public long get_submit_time() {
        return submit_time;
    }

    private void set_submit_time (long submit_time) {
        this.submit_time = submit_time;
    }

    public String get_submit_id() {
        return submit_id;
    }

    private void set_submit_id (String submit_id) {
        this.submit_id = submit_id;
    }

    public int get_opcode() {
        return opcode;
    }

    private void set_opcode (int opcode) {
        this.opcode = opcode;
    }

    public int get_stage() {
        return stage;
    }

    private void set_stage (int stage) {
        this.stage = stage;
    }

    public int get_rescode() {
        return rescode;
    }

    private void set_rescode (int rescode) {
        this.rescode = rescode;
    }

    public String get_results() {
        return results;
    }

    private void set_results (String results) {
        this.results = results;
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

		throw new IllegalArgumentException("LogEntry.set_details: Incorrect type of marshal writer");
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
		String str = "LogEntry\n"
			+ "\tid: " + ((id == null) ? ("null") : (id.toHexString())) + "\n"
			+ "\tlog_time: " + log_time + "\n"
			+ "\tevent_id: " + event_id + "\n"
			+ "\tsched_time: " + sched_time + "\n"
			+ "\tsubmit_time: " + submit_time + "\n"
			+ "\tsubmit_id: " + submit_id + "\n"
			+ "\topcode: " + opcode + "\n"
			+ "\tstage: " + stage + "\n"
			+ "\tdetails: " + get_details_description() + "\n"
			+ "\trescode: " + rescode + "\n"
			+ "\tresults: " + results;
		return str;
	}




	/**
	 * get_record_key - Get the record key for this log entry.
	 */
	public RecordKey get_record_key () {
		return new RecordKey(id);
	}




	/**
	 * set_record_key - Set the record key for this log entry.
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
	 * submit_log_entry - Submit a log entry.
	 * @param key = Record key associated with this task. Can be null to assign a new one.
	 * @param log_time = Time of this log entry, in milliseconds
	 *                   since the epoch. Must be positive.
	 * @param event_id = Event associated with this task, or "" if none. Cannot be null.
	 * @param sched_time = Time at which task should execute, in milliseconds
	 *                     since the epoch. Must be positive.
	 * @param submit_time = Time at which the task is submitted, in milliseconds
	 *                      since the epoch. Must be positive.
	 * @param submit_id = Person or entity submitting this task. Cannot be empty or null.
	 * @param opcode = Operation code used to dispatch the task.
	 * @param stage = Stage number, user-defined, effectively an extension of the opcode.
	 * @param details = Further details of this task. Can be null if there are none.
	 * @param rescode = Result code.
	 * @param results = Further results of this task, or "" if none. Cannot be null.
	 * @return
	 * Returns the new entry.
	 */
	public static LogEntry submit_log_entry (RecordKey key, long log_time, String event_id,
			long sched_time, long submit_time, String submit_id, int opcode, int stage,
			MarshalWriter details, int rescode, String results) {

		// Check conditions

		if (!( log_time > 0L
			&& event_id != null
			&& sched_time > 0L
			&& submit_time > 0L
			&& submit_id != null && submit_id.length() > 0
			&& results != null )) {
			throw new IllegalArgumentException("LogEntry.submit_log_entry: Invalid log parameters");
		}

		// Construct the log entry object

		LogEntry lentry = new LogEntry();
		lentry.set_record_key (key);
		lentry.set_log_time (log_time);
		lentry.set_event_id (event_id);
		lentry.set_sched_time (sched_time);
		lentry.set_submit_time (submit_time);
		lentry.set_submit_id (submit_id);
		lentry.set_opcode (opcode);
		lentry.set_stage (stage);
		lentry.set_details (details);
		lentry.set_rescode (rescode);
		lentry.set_results (results);

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(lentry);
		
		return lentry;
	}




	/**
	 * submit_log_entry - Submit a log entry.
	 * @param ptask = Pending task record associated with this task. Cannot be null.
	 * @param log_time = Time of this log entry, in milliseconds
	 *                     since the epoch. Must be positive.
	 * @param rescode = Result code.
	 * @param results = Further results of this task, or "" if none. Cannot be null.
	 * @return
	 * Other log parameters are copied from ptask.
	 * Returns the new entry.
	 */
	public static LogEntry submit_log_entry (PendingTask ptask, long log_time, int rescode, String results) {

		// Check conditions

		if (!( ptask != null
			&& log_time > 0L
			&& results != null )) {
			throw new IllegalArgumentException("LogEntry.submit_log_entry: Invalid log parameters");
		}

		// Submit the log entry

		LogEntry lentry = submit_log_entry (
			ptask.get_record_key(),
			log_time,
			ptask.get_event_id(),
			ptask.get_sched_time(),
			ptask.get_submit_time(),
			ptask.get_submit_id(),
			ptask.get_opcode(),
			ptask.get_stage(),
			ptask.get_details_as_payload(),
			rescode,
			results);
		
		return lentry;
	}




	/**
	 * store_log_entry - Store a log entry into the database.
	 * This is primarily for restoring from backup.
	 */
	public static LogEntry store_log_entry (LogEntry lentry) {

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(lentry);
		
		return lentry;
	}




	/**
	 * get_log_entry_for_key - Get the log entry with the given key.
	 * @param key = Record key. Cannot be null or empty.
	 * Returns the log entry, or null if not found.
	 */
	public static LogEntry get_log_entry_for_key (RecordKey key) {

		if (!( key != null && key.getId() != null )) {
			throw new IllegalArgumentException("LogEntry.get_log_entry_for_key: Missing or empty record key");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select id == key

		Query<LogEntry> query = datastore.createQuery(LogEntry.class)
											.filter("id ==", key.getId());

		// Run the query

		LogEntry lentry = query.get();

		return lentry;
	}




	/**
	 * get_log_entry_range - Get a range of log entries, reverse-sorted by log time.
	 * @param log_time_lo = Minimum log time, in milliseconds since the epoch.
	 *                      Can be 0L for no minimum.
	 * @param log_time_hi = Maximum log time, in milliseconds since the epoch.
	 *                      Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 */
	public static List<LogEntry> get_log_entry_range (long log_time_lo, long log_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<LogEntry> query = datastore.createQuery(LogEntry.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with log_time >= log_time_lo

		if (log_time_lo > 0L) {
			query = query.filter("log_time >=", new Long(log_time_lo));
		}

		// Select entries with log_time <= log_time_hi

		if (log_time_hi > 0L) {
			query = query.filter("log_time <=", new Long(log_time_hi));
		}

		// Sort by log_time in descending order (most recent first)

		query = query.order("-log_time");

		// Run the query

		List<LogEntry> entries = query.asList();

		return entries;
	}




	/**
	 * fetch_log_entry_range - Iterate a range of log entries, reverse-sorted by log time.
	 * @param log_time_lo = Minimum log time, in milliseconds since the epoch.
	 *                      Can be 0L for no minimum.
	 * @param log_time_hi = Maximum log time, in milliseconds since the epoch.
	 *                      Can be 0L for no maximum.
	 * @param event_id = Event id. Can be null to return entries for all events.
	 */
	public static RecordIterator<LogEntry> fetch_log_entry_range (long log_time_lo, long log_time_hi, String event_id) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<LogEntry> query = datastore.createQuery(LogEntry.class);

		// Select by event_id

		if (event_id != null) {
			query = query.filter("event_id ==", event_id);
		}

		// Select entries with log_time >= log_time_lo

		if (log_time_lo > 0L) {
			query = query.filter("log_time >=", new Long(log_time_lo));
		}

		// Select entries with log_time <= log_time_hi

		if (log_time_hi > 0L) {
			query = query.filter("log_time <=", new Long(log_time_hi));
		}

		// Sort by log_time in descending order (most recent first)

		query = query.order("-log_time");

		// Run the query

		MorphiaIterator<LogEntry, LogEntry> morphia_iterator = query.fetch();

		return new RecordIterator<LogEntry>(morphia_iterator);
	}




	/**
	 * delete_log_entry - Delete a log entry.
	 * @param entry = Existing log entry to delete.
	 * @return
	 */
	public static void delete_log_entry (LogEntry entry) {

		// Check conditions

		if (!( entry != null && entry.get_id() != null )) {
			throw new IllegalArgumentException("LogEntry.delete_log_entry: Invalid parameters");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Run the delete

		datastore.delete(entry);
		
		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 9001;

	private static final String M_VERSION_NAME = "LogEntry";

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		String sid = id.toHexString();
		writer.marshalString      ("id"         , sid        );
		writer.marshalLong        ("log_time"   , log_time   );
		writer.marshalString      ("event_id"   , event_id   );
		writer.marshalLong        ("sched_time" , sched_time );
		writer.marshalLong        ("submit_time", submit_time);
		writer.marshalString      ("submit_id"  , submit_id  );
		writer.marshalInt         ("opcode"     , opcode     );
		writer.marshalInt         ("stage"      , stage      );
		writer.marshalJsonString  ("details"    , details    );
//		writer.marshalLongArray   ("details_l"  , details_l  );
//		writer.marshalDoubleArray ("details_d"  , details_d  );
//		writer.marshalStringArray ("details_s"  , details_s  );
		writer.marshalInt         ("rescode"    , rescode    );
		writer.marshalString      ("results"    , results    );
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		String sid;
		sid         = reader.unmarshalString      ("id"         );
		log_time    = reader.unmarshalLong        ("log_time"   );
		event_id    = reader.unmarshalString      ("event_id"   );
		sched_time  = reader.unmarshalLong        ("sched_time" );
		submit_time = reader.unmarshalLong        ("submit_time");
		submit_id   = reader.unmarshalString      ("submit_id"  );
		opcode      = reader.unmarshalInt         ("opcode"     );
		stage       = reader.unmarshalInt         ("stage"      );
		details     = reader.unmarshalJsonString  ("details"    );
//		details_l   = reader.unmarshalLongArray   ("details_l"  );
//		details_d   = reader.unmarshalDoubleArray ("details_d"  );
//		details_s   = reader.unmarshalStringArray ("details_s"  );
		rescode     = reader.unmarshalInt         ("rescode"    );
		results     = reader.unmarshalString      ("results"    );
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

	public LogEntry unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}


}
