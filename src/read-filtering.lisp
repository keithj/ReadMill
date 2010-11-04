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

(defun quality-filter-bam (cmd input-file output-file start end
                           quality-threshold)
  "Filters BAM file by read quality, removing any reads that have a
base quality less than QUALITY-THRESHOLD between base positions START
and END. Updates the BAM header in the output file with a new PG
header record."
  (with-bam (in (header num-refs ref-meta) input-file)
    (let* ((hd (make-sam-header header))
           (pg (pg-record "readmill-quality-filter"
                          :program-name *readmill-name*
                          :program-version *readmill-version*
                          :previous-program (previous-program hd)
                          :command-line (format nil "~{~a~^ ~}" cmd))))
      (with-bam (out ((with-output-to-string (s)
                        (write-sam-header
                         (add-pg-record hd pg) s)) num-refs ref-meta) output-file
                         :direction :output
                         :if-does-not-exist :create :if-exists :overwrite)
        (let* ((start (or start 0))
               (qfilter (make-counting-predicate
                         (make-quality-p quality-threshold
                                         :start start :end end)))
               (in (discarding-if qfilter in)))
          (loop
             while (has-more-p in)
             do (consume out (next in)))
          (rest (funcall qfilter nil)))))))

(defun subseq-filter-bam (cmd input-file output-file start end queries)
  "Filters BAM file by read subsequence, removing any reads that have
any of the strings in list QUERIES between base positions START and
END. Updates the BAM header in the output file with a new PG header
record."
  (with-bam (in (header num-refs ref-meta) input-file)
    (let* ((hd (make-sam-header header))
           (pg (pg-record "readmill-subseq-filter"
                          :program-name *readmill-name*
                          :program-version *readmill-version*
                          :previous-program (previous-program hd)
                          :command-line (format nil "~{~a~^ ~}" cmd))))
      (with-bam (out ((with-output-to-string (s)
                        (write-sam-header
                         (add-pg-record hd pg) s)) num-refs ref-meta) output-file
                         :direction :output
                         :if-does-not-exist :create :if-exists :overwrite)
        (let* ((start (or start 0))
               (fns (mapcar (lambda (query)
                              (make-counting-predicate
                               (make-subseq-p query :start start :end end)))
                            queries))
               (in (discarding-if (apply #'any fns) in)))
          (loop
             while (has-more-p in)
             do (consume out (next in)))
          (loop
             for fn in fns
             collect (rest (multiple-value-list (funcall fn nil)))))))))

;; Experimental multi-threaded version
(defun subseq-pfilter-bam (cmd input-file output-file start end queries)
  (with-bam (in (header num-refs ref-meta) input-file)
    (let* ((hd (make-sam-header header))
           (pg (pg-record "readmill-subseq-filter"
                          :program-name *software-name*
                          :program-version *software-version*
                          :previous-program (previous-program hd)
                          :command-line (format nil "~{~a~^ ~}" cmd))))
      (with-bam (out ((with-output-to-string (s)
                        (write-sam-header
                         (add-pg-record hd pg) s)) num-refs ref-meta)
                     output-file :direction :output
                     :if-does-not-exist :create :if-exists :overwrite)
        (let* ((start (or start 0))
               (fns (mapcar (lambda (query)
                              (make-subseq-p query :start start :end end))
                            queries)))
          (batch-multi-filter in out fns :threads 2 :batch-size 10000))))))
