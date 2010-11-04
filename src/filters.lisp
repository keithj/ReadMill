;;;
;;; Copyright (C) 2010 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of readmill.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :uk.ac.sanger.readmill)

(defmacro with-bam ((var (&optional header num-refs ref-meta) filespec
                         &rest args &key (compress t) (null-padding 0)
                         &allow-other-keys)
                    &body body)
  "Evaluates BODY with VAR bound to a newly opened BAM iterator on
pathname designator FILESPEC. The direction (:input versus :output) is
determined by the stream-opening arguments in ARGS.

On reading, HEADER, NUM-REFS and REF-META will be automatically bound
to the BAM file metadata i.e. the metadata are read automatically and
the iterator positioned before the first alignment record.

On writing, HEADER, NUM-REFS and REF-META should be bound to
appropriate values for the BAM file metadata, which will be
automatically written to the underlying stream.

The COMPRESS and NULL-PADDING keyword arguments are only applicable on
writing where they control whether the header block should be
compressed and whether the header string should be padded with nulls
to allow space for expansion.

For example:

To count all records in a BAM file:

;;; (with-bam (in (header) \"in.bam\")
;;;   (declare (ignore header))
;;;   (loop
;;;      while (has-more-p in)
;;;      count (next in)))

To copy only records with a mapping quality of >= 30 to another BAM
file:

;;; (with-bam (in (h n r) \"in.bam\")
;;;   (with-bam (out (h n r) \"out.bam\" :direction :output)
;;;     (let ((q30 (discarding-if (lambda (x)
;;;                                 (< (mapping-quality x) 30)) in)))
;;;       (loop
;;;          while (has-more-p q30)
;;;          do (consume out (next q30))))))"
  (with-gensyms (bgzf)
    (let ((default-header (with-output-to-string (s)
                            (write-sam-header
                             `((:HD (:VN . ,*sam-version*))) s))))
      `(with-bgzf (,bgzf ,filespec ,@(remove-key-values
                                      '(:compress :null-padding) args))
         ,@(if (search '(:direction :output) args)
               `((write-bam-meta ,bgzf
                                 ,(or header default-header)
                                 ,(or num-refs 0) ,ref-meta
                                 :compress ,compress
                                 :null-padding ,null-padding)
                 (let ((,var (make-bam-output ,bgzf)))
                   ,@body))
               `((multiple-value-bind (,@(when header `(,header))
                                       ,@(when num-refs `(,num-refs))
                                       ,@(when ref-meta `(,ref-meta)))
                     (read-bam-meta ,bgzf)
                   ,@(cond (ref-meta
                            `((declare (ignorable ,header ,num-refs))))
                           (num-refs
                            `((declare (ignorable ,header)))))
                   (let ((,var (make-bam-input ,bgzf)))
                     ,@body))))))))

(defun make-counting-predicate (predicate)
  "Returns a copy of PREDICATE that counts the number of alignments
filtered and the number passed. PREDICATE should accept a BAM
record. The returned function returns three values; the return value
of PREDICATE, the number of times the function has been called with a
BAM record argument and the number of times that it has returned T."
  (let ((called 0)
        (passed 0))
    (declare (type fixnum called passed))
    (lambda (&optional x)
      (declare (optimize (speed 3) (safety 0)))
      (declare (type function predicate))
      (cond ((null x)
             (values x called passed))
            ((funcall predicate x)
             (incf called)
             (values x called passed))
            (t
             (incf called)
             (incf passed)
             (values nil called passed))))))

(defun make-bam-input (bam)
  "Returns a generator function that returns BAM alignment records for
BAM stream BAM. The standard generator interface functions NEXT and
HAS-MORE-P may be used in operations on the returned generator.

For example:

To count all records in a BAM file:

;;; (with-bgzf (bam \"in.bam\")
;;;   (read-bam-meta bam)
;;;   (let ((in (make-bam-input bam)))
;;;     (loop
;;;        while (has-more-p in)
;;;        count (next in))))

To copy only records with a mapping quality of >= 30 to another BAM
file:

;;; (with-bgzf (bam-in \"in.bam\")
;;;   (with-bgzf (bam-out \"out.bam\" :direction :output)
;;;     (multiple-value-bind (header num-refs ref-meta)
;;;         (read-bam-meta bam-in)
;;;       (write-bam-meta bam-out header num-refs ref-meta))
;;;     (let ((fn (discarding-if (lambda (x)
;;;                                (< (mapping-quality x) 30))
;;;                              (make-bam-input bam-in))))
;;;       (loop
;;;          while (has-more-p fn)
;;;          do (write-alignment bam-out (next fn))))))"
  (let ((current (read-alignment bam)))
    (defgenerator
        (more (not (null current)))
        (next (prog1
                  current
                (setf current (read-alignment bam)))))))

(defun make-bam-output (bam)
  "Returns a consumer function that accepts an argument of a BAM
record and writes it to BAM output stream BAM. The standard consumer
interface function CONSUME may be used in operations on the returned
consumer."
  (lambda (aln)
    (write-alignment bam aln)))

(defun batch-filter (in out predicate &key (threads 1) (batch-size 1000))
  "Obtains values from generator function IN and passes any for which
function PREDICATE does not return T, to consumer function
OUT (i.e. acts like CL:DELETE-IF). Calls to PREDICATE may be run in
batches under separate threads, each processing BATCH-SIZE records.
The return values are a list of the filter predicates created and a
list of the counts of values tested and passed."
  (declare (optimize (speed 3)))
  (declare (type fixnum threads batch-size))
  (flet ((collect-batches (m n)
           (let ((batches ()))
             (dotimes (x m (nreverse batches))
               (push (collect in n) batches))))
         (process-batch (fn batch)
           (declare (type function fn)
                    (type list batch))
           (eager-future:pexec (delete-if fn batch))))
    (let ((fns (mapcar #'make-counting-predicate (loop
                                                    repeat threads
                                                    collect predicate))))
      (loop
         for batches of-type list = (collect-batches threads batch-size)
         until (every #'null batches)
         do (let ((futures (mapcar #'process-batch fns batches)))
              (dolist (result (mapcar #'eager-future:yield futures))
                (dolist (x result)
                  (consume out x)))))
      (mapcar (lambda (fn)
                (declare (type function fn))
                (rest (multiple-value-list (funcall fn)))) fns))))

;; Experimental version that uses a list of predicates that are
;; wrapped in in counters, which are in turn composed into one
;; function per thread
(defun batch-multi-filter (in out predicates &key (threads 1) (batch-size 1000))
  (declare (optimize (speed 3)))
  (declare (type fixnum threads batch-size))
  (flet ((collect-batches (m n)
           (let ((batches ()))
             (dotimes (x m (nreverse batches))
               (push (collect in n) batches))))
         (process-batch (fn batch)
           (declare (type function fn)
                    (type list batch))
           (eager-future:pexec (delete-if fn batch)))
         (collect-counts (counters)
           (mapcar (lambda (fn)
                     (declare (type function fn))
                     (rest (multiple-value-list (funcall fn)))) counters)))
    (let* ((counter-sets (loop
                            repeat threads
                            collect (mapcar #'make-counting-predicate
                                            predicates)))
           (fns (mapcar (lambda (counters)
                          (apply #'any counters)) counter-sets)))
      (loop
         for batches of-type list = (collect-batches threads batch-size)
         until (every #'null batches)
         do (let ((futures (mapcar #'process-batch fns batches)))
              (dolist (passed (mapcar #'eager-future:yield futures))
                (dolist (obj passed)
                  (consume out obj)))))
      (mapcar #'collect-counts counter-sets))))

(defun make-rg-p (read-group)
  "Returns a read group filtering predicate that returns T for BAM
alignments of READ-GROUP."
  (lambda (aln)
    (declare (optimize (speed 3)))
    (declare (type simple-base-string read-group))
    (let ((tag (assocdr :rg (alignment-tag-values aln))))
      (when tag
        (locally (declare (type simple-base-string tag))
          (string/= read-group tag))))))

(defun make-polyx-p (char length &key (start 0) end)
  "Returns a homopolymer filtering predicate that returns T for BAM
alignments containing runs of CHAR or LENGTH or longer, between START
and END of the alignment unclipped sequence string."
  (make-subseq-p (make-string length
                              :element-type 'base-char
                              :initial-element char) :start start :end end))

(defun make-subseq-p (str &key (start 0) end)
  "Returns a subsequence matching predicate that returns T for BAM
alignments containing subsequence STR between START and END of the
alignment unclipped sequence string."
  (let ((str (string-upcase str)))
    (lambda (aln)
      (declare (optimize (speed 3)))
      (let ((seq (seq-string aln)))
        (declare (type simple-base-string seq))
        (search str seq :start2 start :end2 end :test #'string=)))))

(defun make-quality-p (quality-threshold &rest args
                       &key (start 0) end (test #'<=) (count 1))
  "Returns a quality matching predicate that returns T for BAM
alignments having COUNT or more base qualities passing TEST relative
to QUALITY-THRESHOLD between READ-START and READ-END. Test defaults to
<= and COUNT defaults to 1, i.e. the default behaviour is to return T
where at least 1 quality value is <= QUALITY-THRESHOLD between START
and END."
  (declare (ignorable start end test))
  (check-arguments (<= 0 quality-threshold 100) (quality-threshold)
                   "Quality threshold was not in the range 0-100")
  (check-arguments (and (integerp count) (plusp count)) (count)
                   "expected a positive integer")
  (lambda (aln)
    (declare (optimize (speed 3)))
    (let ((qual-str (quality-string aln)))
      (declare (type simple-base-string qual-str)
               (type quality-score quality-threshold)
               (type fixnum count))
      (>= count (apply #'count quality-threshold
                       (map-into (make-array (length qual-str)
                                             :element-type 'quality-score
                                             :initial-element 0)
                                 #'decode-phred-quality qual-str)
                       (remove-key-values '(:count) args))))))

(defun maybe-filter-rg (in header &optional read-group)
  (if read-group
      (let ((rg (find read-group (header-records (make-sam-header header) :rg)
                      :key (lambda (record)
                             (header-value record :id))
                      :test #'string=)))
        (check-arguments rg (read-group) "this read-group is not present")
        (values
         (discarding-if (lambda (aln)
                          (string/= read-group
                                    (assocdr :rg (alignment-tag-values aln))))
                        in)
         (format nil "read group ~a ~a~@[ (~a)~]" read-group
                 (header-value rg :sm) (header-value rg :pu))))
      (values in "all read groups")))
