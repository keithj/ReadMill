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

;; (defmacro with-bam-input ((var (header &optional num-refs ref-meta)
;;                                filespec &rest args) &body body)
;;   (with-gensyms (bgzf)
;;     `(with-bgzf (,bgzf ,filespec ,@args)
;;        (multiple-value-bind (,header ,@(when num-refs `(,num-refs))
;;                                      ,@(when ref-meta `(,ref-meta)))
;;            (read-bam-meta ,bgzf)
;;          (let ((,var (make-bam-input ,bgzf)))
;;            ,@body)))))

;; (defmacro with-bam-output ((var (header num-refs ref-meta &key (compress t)
;;                                         (null-padding 0)) filespec &rest args)
;;                            &body body)
;;   (with-gensyms (bgzf)
;;     `(with-bgzf (,bgzf ,filespec ,@args)
;;        (write-bam-meta ,bgzf ,header ,num-refs ,ref-meta :compress ,compress
;;                        :null-padding ,null-padding)
;;        (let ((,var (make-bam-output ,bgzf)))
;;          ,@body))))

(defmacro with-bam ((var (header &optional num-refs ref-meta) filespec
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
    `(with-bgzf (,bgzf ,filespec ,@(remove-key-values
                                    '(:compress :null-padding) args))
       ,@(if (search '(:direction :output) args)
             `((write-bam-meta ,bgzf ,header ,num-refs ,ref-meta
                               :compress ,compress
                               :null-padding ,null-padding)
               (let ((,var (make-bam-output ,bgzf)))
                 ,@body))
             `((multiple-value-bind (,header ,@(when num-refs `(,num-refs))
                                             ,@(when ref-meta `(,ref-meta)))
                   (read-bam-meta ,bgzf)
                 (let ((,var (make-bam-input ,bgzf)))
                   ,@body)))))))

(defun make-counting-predicate (predicate)
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
    "Obtains values from generator function IN and passes any for
which function PREDICATE does not return T, to consumer function
OUT (i.e. acts like CL:DELETE-IF). Calls to PREDICATE may be run in
batches under separate threads, each processing BATCH-SIZE records.

The return values are a list of the filter predicates created and a
list of the counts of values tested and passed."
    (declare (optimize (speed 3)))
    (declare (type fixnum threads batch-size))
    (let ((fns (mapcar #'make-counting-predicate (loop
                                                    repeat threads
                                                    collect predicate))))
      (loop
         for batches of-type list = (collect-mxn in threads batch-size)
         until (every #'null batches)
         do (let ((futures (mapcar (lambda (fn batch)
                                     (declare (type function fn)
                                              (type list batch))
                                     (eager-future:pexec
                                      (delete-if fn batch))) fns batches)))
              (dolist (result (mapcar #'eager-future:yield futures))
                (dolist (x result)
                  (consume out x)))))
      (values fns (mapcar (lambda (fn)
                            (declare (type function fn))
                            (multiple-value-list (funcall fn))) fns))))

(defun collect-mxn (in m n)
  "Collects M lists of up to N records each from BAM producer function
IN."
  (declare (optimize (speed 3)))
  (declare (type fixnum m n))
  (let ((batches ()))
    (dotimes (x m (nreverse batches))
      (push (collect in n) batches))))

