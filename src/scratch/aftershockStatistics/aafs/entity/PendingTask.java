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

import scratch.aftershockStatistics.aafs.MongoDBUtil;
import scratch.aftershockStatistics.aafs.RecordKey;

import scratch.aftershockStatistics.MarshalImpArray;
import scratch.aftershockStatistics.MarshalReader;
import scratch.aftershockStatistics.MarshalWriter;




/**
 * Holds a pending AAFS task in the MongoDB database.
 * Author: Michael Barall 03/15/2018.
 *
 * The collection "tasks" holds the queue of pending tasks.
 */
@Entity(value = "tasks", noClassnameStored = true)
@Indexes({
	@Index(fields = {@Field("event_id")}, options = @IndexOptions(name = "eventid")),
	@Index(fields = {@Field("exec_time")}, options = @IndexOptions(name = "extime"))
})
public class PendingTask implements java.io.Serializable {

	//----- Envelope information -----

	// Globally unique identifier for this task.
	// This is the MongoDB identifier.
	// Note that ObjectId implements java.io.Serializable.

    @Id
    private ObjectId id;

	// Time that this task is scheduled to execute, in milliseconds since the epoch.
	// The collection is indexed on this field, so the task with the earliest
	// time can be obtained with a database query efficiently.
	// When beginning to execute a task, this field is set to zero, which guarantees
	// it will be seen again if the task is interrupted and restarted.

	//@Indexed(options = @IndexOptions(name = "extime"))
	private long exec_time;

	//----- Task information -----

	// Event ID that this task pertains to.
	// Tasks not referring to an event should put an empty string here (not null).
	// This can be used to find all queued tasks pertaining to an event.

	//@Indexed(options = @IndexOptions(name = "eventid"))
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
	// Any additional information neeed is stored as marshaled data.
	// Each array should have at least one element.

	@Embedded
	private long[] details_l;
	@Embedded
	private double[] details_d;
	@Embedded
	private String[] details_s;




	//----- Getters and setters -----

    private ObjectId get_id() {
        return id;
    }

    private void set_id (ObjectId id) {
        this.id = id;
    }

    private long get_exec_time() {
        return exec_time;
    }

    private void set_exec_time (long exec_time) {
        this.exec_time = exec_time;
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




	/**
	 * get_details - Get a reader for the details.
	 */
    public MarshalReader get_details() {
        return new MarshalImpArray (details_l, details_d, details_s);
    }


	/**
	 * set_details - Set details from the marshaled data, can be null for none.
	 */
    private void set_details (MarshalWriter writer) {
		if (writer == null) {
			details_l = new long[1];
			details_l[0] = 0L;
			details_d = new double[1];
			details_d[0] = 0.0;
			details_s = new String[1];
			details_s[0] = "";
			return;
		}

		if (!( writer instanceof MarshalImpArray )) {
			throw new IllegalArgumentException("PendingTask.set_details: Incorrect type of marshal writer");
		}

		MarshalImpArray w = (MarshalImpArray)writer;
		details_l = w.get_long_store();
		details_d = w.get_double_store();
		details_s = w.get_string_store();
		return;
    }


	/**
	 * begin_details - Get a writer to use for marhaling details.
	 */
    public static MarshalWriter begin_details() {
        return new MarshalImpArray ();
    }


	/**
	 * get_details_as_writer - Get a writer containing the details.
	 */
    MarshalWriter get_details_as_writer() {
        return new MarshalImpArray (details_l, details_d, details_s);
    }


	///**
	// * These getters are for use by other package members for direct copy of the details.
	// */
    //long[] get_details_l() {
    //    return details_l;
    //}
    //double[] get_details_d() {
    //    return details_d;
    //}
    //String[] get_details_s() {
    //    return details_s;
    //}




	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "PendingTask\n"
			+ "\tid: " + ((id == null) ? ("null") : (id.toHexString())) + "\n"
			+ "\texec_time: " + exec_time + "\n"
			+ "\tevent_id: " + event_id + "\n"
			+ "\tsched_time: " + sched_time + "\n"
			+ "\tsubmit_time: " + submit_time + "\n"
			+ "\tsubmit_id: " + submit_id + "\n"
			+ "\topcode: " + opcode + "\n"
			+ "\tstage: " + stage + "\n"
			+ "\tdetails: " + "llen = " + details_l.length + ", dlen = " + details_d.length + ", slen = " + details_s.length;
		return str;
	}




	/**
	 * get_record_key - Get the record key for this pending task.
	 * Each pending task is assigned a unique record key when the task is submitted.
	 * The key remains the same when the task is activated or when a new stage begins.
	 */
	public RecordKey get_record_key () {
		return new RecordKey(id);
	}




	/**
	 * is_restarted - Return true if this task has been restarted.
	 * Restarted means that a previous attempt was made to execute the task,
	 * but it failed before deleting the task from the queue.
	 */
	public boolean is_restarted () {
		return exec_time == 0L;
	}




	/**
	 * get_apparent_time - Return the apparent execution time of this command.
	 * Note: This is primarily for accelerated testing.
	 */
	public long get_apparent_time () {
		return (exec_time == 0L) ? sched_time : exec_time;
	}




	/**
	 * submit_task - Submit a task.
	 * @param event_id = Event associated with this task, or "" if none. Cannot be null.
	 * @param sched_time = Time at which task should execute, in milliseconds
	 *                     since the epoch. Must be positive.
	 * @param submit_time = Time at which the task is submitted, in milliseconds
	 *                      since the epoch. Must be positive.
	 * @param submit_id = Person or entity submitting this task. Cannot be empty or null.
	 * @param opcode = Operation code used to dispatch the task.
	 * @param stage = Stage number, user-defined, effectively an extension of the opcode.
	 * @param details = Further details of this task. Can be null if there are none.
	 * @return
	 * Returns the new entry.
	 */
	public static PendingTask submit_task (String event_id, long sched_time, long submit_time,
								String submit_id, int opcode, int stage, MarshalWriter details) {

		// Check conditions

		if (!( event_id != null
			&& sched_time > 0L
			&& submit_time > 0L
			&& submit_id != null && submit_id.length() > 0 )) {
			throw new IllegalArgumentException("PendingTask.submit_task: Invalid task parameters");
		}

		// Construct the pending task object

		PendingTask ptask = new PendingTask();
		ptask.set_id (null);
		ptask.set_exec_time (sched_time);
		ptask.set_event_id (event_id);
		ptask.set_sched_time (sched_time);
		ptask.set_submit_time (submit_time);
		ptask.set_submit_id (submit_id);
		ptask.set_opcode (opcode);
		ptask.set_stage (stage);
		ptask.set_details (details);

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(ptask);
		
		return ptask;
	}




	/**
	 * get_all_tasks_unsorted - Get a list of all pending tasks, without sorting.
	 * This is primarily for testing and monitoring.
	 */
	public static List<PendingTask> get_all_tasks_unsorted () {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select documents of type PendingTask

		Query<PendingTask> query = datastore.createQuery(PendingTask.class);

		// Run the query

		List<PendingTask> tasks = query.asList();

		return tasks;
	}




	/**
	 * get_all_tasks - Get a list of all pending tasks, sorted by execution time.
	 * This is primarily for testing and monitoring.
	 */
	public static List<PendingTask> get_all_tasks () {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Sort by exec_time

		Query<PendingTask> query = datastore.createQuery(PendingTask.class)
											.order("exec_time");

		// Run the query

		List<PendingTask> tasks = query.asList();

		return tasks;
	}




	/**
	 * get_first_task - Get the first task, that is, the task with smallest execution time.
	 * This is primarily for testing and monitoring.
	 */
	public static PendingTask get_first_task () {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Sort by exec_time

		Query<PendingTask> query = datastore.createQuery(PendingTask.class)
											.order("exec_time");

		// Run the query

		PendingTask task = query.get();

		return task;
	}




	/**
	 * get_first_ready_task - Get the first ready task, according to execution time.
	 * @param cutoff_time = Cutoff time, in milliseconds since the epoch.
	 * Only tasks with exec_time <= cutoff_time are considered.
	 * Return is null if there are no such tasks.
	 * This is primarily for testing and monitoring.
	 */
	public static PendingTask get_first_ready_task (long cutoff_time) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select exec_time <= cutoff_time, sort by exec_time

		Query<PendingTask> query = datastore.createQuery(PendingTask.class)
											.filter("exec_time <=", new Long(cutoff_time))
											.order("exec_time");

		// Run the query

		PendingTask task = query.get();

		return task;
	}




	/**
	 * activate_first_ready_task - Get and activate the first ready task, according to execution time.
	 * @param cutoff_time = Cutoff time, in milliseconds since the epoch.
	 * Only tasks with exec_time <= cutoff_time are considered.
	 * Return is null if there are no such tasks.
	 * The task is marked active by setting exec_time = 0 in the database.
	 */
	public static PendingTask activate_first_ready_task (long cutoff_time) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select exec_time <= cutoff_time, sort by exec_time

		Query<PendingTask> query = datastore.createQuery(PendingTask.class)
											.filter("exec_time <=", new Long(cutoff_time))
											.order("exec_time");

		// Construct the update operation: Set exec_time to 0L

		UpdateOperations<PendingTask> update_op
				= datastore.createUpdateOperations(PendingTask.class)
							.set("exec_time", new Long(0L));

		// Construct the find and modify options: Return the original document value

		FindAndModifyOptions modify_opt = (new FindAndModifyOptions())
											.returnNew(false);

		// Run the query

		PendingTask task = datastore.findAndModify(query, update_op, modify_opt);

		return task;
	}




	/**
	 * stage_task - Begin a new stage of a task.
	 * @param ptask = Existing pending task to stage.
	 * @param exec_time = Time at which task should execute, in milliseconds
	 *                    since the epoch. Must be positive.
	 * @param stage = Stage number, user-defined, effectively an extension of the opcode.
	 * @return
	 */
	public static void stage_task (PendingTask ptask, long exec_time, int stage) {

		// Check conditions

		if (!( ptask != null && ptask.get_id() != null
			&& exec_time > 0L )) {
			throw new IllegalArgumentException("PendingTask.stage_task: Invalid task parameters");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the update operations: Set exec_time and stage

		UpdateOperations<PendingTask> update_op
				= datastore.createUpdateOperations(PendingTask.class)
							.set("exec_time", new Long(exec_time))
							.set("stage", new Integer(stage));

		// Run the update

		datastore.update(ptask, update_op);
		
		return;
	}




	/**
	 * delete_task - Delete a task.
	 * @param ptask = Existing pending task to delete.
	 * @return
	 */
	public static void delete_task (PendingTask ptask) {

		// Check conditions

		if (!( ptask != null && ptask.get_id() != null )) {
			throw new IllegalArgumentException("PendingTask.delete_task: Invalid task parameters");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Run the delete

		datastore.delete(ptask);
		
		return;
	}




}
