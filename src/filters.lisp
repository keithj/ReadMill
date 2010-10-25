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

(defun make-filter-predicate (predicate)
  "Returns a copy of PREDICATE that counts the number of alignments
filtered and the number passed. PREDICATE should accept a BAM
record. The returned function returns three values; the return value
of PREDICATE, the number of times the function has been called with a
BAM record argument and the number of times that it has returned T."
  (let ((filtered 0)
        (passed 0))
    (declare (type fixnum filtered passed))
    (lambda (&optional x)
      (declare (optimize (speed 3) (safety 0)))
      (declare (type function predicate))
      (cond ((null x)
             (values nil filtered passed))
            ((funcall predicate x)
             (incf filtered)
             (incf passed)
             (values t filtered passed))
            (t
             (incf filtered)
             (values nil filtered passed))))))

(defun collect-mxn (bam m n)
  "Collects M lists of up to N records from BAM stream."
  (declare (optimize (speed 3)))
  (declare (type fixnum m n))
  (loop
     repeat m
     collect (loop
                repeat n
                for aln = (read-alignment bam)
                while aln
                collect aln)))

(defun filter-bam (bam1 bam2 predicate &key (threads 1) (batch-size 1000))
  "Reads BAM records from BAM1 and writes any for which function
PREDICATE returns T, to BAM2. Calls to PREDICATE may be run in batches
under separate threads, each processing BATCH-SIZE records.

The return values are a list of the filter predicates created and a
list of the counts of BAM records tested and passed."
  (declare (optimize (speed 3)))
  (declare (type fixnum threads batch-size))
  (let ((fns (mapcar #'make-filter-predicate
                     (loop
                        repeat threads
                        collect predicate))))
    (loop
       for batches of-type list = (collect-mxn bam1 threads batch-size)
       until (every #'null batches)
       do (let ((futures (mapcar (lambda (fn batch)
                                   (declare (type function fn)
                                            (type list batch))
                                   (eager-future:pexec
                                     (delete-if fn batch))) fns batches)))
            (dolist (result (mapcar #'eager-future:yield futures))
              (dolist (aln result)
                (write-alignment bam2 aln)))))
    (values fns (mapcar (lambda (fn)
                          (declare (type function fn))
                          (multiple-value-list (funcall fn))) fns))))
