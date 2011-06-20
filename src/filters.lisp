;;;
;;; Copyright (c) 2010-2011 Genome Research Ltd. All rights reserved.
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

(defun make-counting-predicate (predicate)
  "Returns a copy of PREDICATE that counts the number of alignments
filtered. PREDICATE should accept a BAM record. The returned function
returns three values; the return value of PREDICATE, the number of
times the function has been called with a BAM record argument and the
number of times that it has returned T."
  (let ((calls 0)
        (true 0))
    (declare (type fixnum calls true))
    (lambda (&optional x)
      (declare (optimize (speed 3) (safety 0)))
      (declare (type function predicate))
      (cond ((null x)
             (values x calls true))
            ((funcall predicate x)
             (incf calls)
             (incf true)
             (values x calls true))
            (t
             (incf calls)
             (values nil calls true))))))

;; Experimental version that uses a list of predicates that are
;; wrapped in in counters, which are in turn composed into one
;; function per thread
(defun batch-filter (in out predicates &key (threads 1) (batch-size 1000))
  "Obtains values from generator function IN and passes any for which
any of functions PREDICATES do not return T, to consumer function
OUT (i.e. acts like CL:DELETE-IF). Calls to PREDICATES may be run in
batches under separate threads, each processing BATCH-SIZE records.
The return value is a list of lists, each containing the counts of
values tested and passed as the first and second elements,
respectively."
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
                     (rest (multiple-value-list (funcall fn)))) counters))
         (aggregate-counts (count-lists)
           (apply #'mapcar (lambda (&rest args)
                             args) count-lists))
         (sum-counts (counts)
           (list (reduce #'+ counts :key #'first)
                 (reduce #'+ counts :key #'second))))
    (let* ((counter-sets (loop
                            repeat threads
                            collect (mapcar #'make-counting-predicate
                                            predicates)))
           (fns (mapcar (lambda (counters)   ; one fn per thread
                          (apply #'any counters)) counter-sets)))
      (loop
         for batches of-type list = (collect-batches threads batch-size)
         until (every #'null batches)
         do (let ((futures (mapcar #'process-batch fns batches)))
              (dolist (passed (mapcar #'eager-future:yield futures))
                (dolist (obj passed)
                  (consume out obj)))))
      (mapcar #'sum-counts
              (aggregate-counts (mapcar #'collect-counts counter-sets))))))

(defun make-rg-p (read-group)
  "Returns a read group filtering predicate that returns T for BAM
alignments of READ-GROUP."
  (lambda (aln)
    (declare (optimize (speed 3)))
    (declare (type simple-string read-group))
    (let ((tag (assocdr :rg (alignment-tag-values aln))))
      (when tag
        (locally (declare (type simple-base-string tag))
          (string= read-group tag))))))

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

(defun make-quality-p (quality-threshold
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
      (>= count (count quality-threshold
                       (map-into (make-array (length qual-str)
                                             :element-type 'quality-score
                                             :initial-element 0)
                                 #'decode-phred-quality qual-str)
                       :start start :end end :test test)))))

(defun pair-consumer (out)
  "Given an alignment consumer function OUT, returns a new consumer
function that only successfully consumes read pairs. i.e. consecutive
first/last fragment alignments in the new template/fragment SAM
notation. Orphan read alignments are dropped."
  (let (aln1 name1)
    (lambda (aln)
      (cond ((and (null aln) (null aln1))
             nil)
            ((and (null aln) aln1)
             (setf aln1 nil))
            (t
             (let ((flag (alignment-flag aln)))
               (cond ((and (null aln1) (multiple-frags-p flag))
                      (setf aln1 aln
                            name1 (read-name aln)))
                     ((string= name1 (read-name aln))
                      (cond ((first-frag-p flag)
                             (consume out aln)
                             (consume out aln1))
                            (t
                             (consume out aln1)
                             (consume out aln)))
                      (setf aln1 nil))
                     (t
                      (setf aln1 aln
                            name1 (read-name aln))))))))))
