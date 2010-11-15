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

(defun filter-bam (cmd input-file output-file filters descriptors
                   &optional json-file)
  (with-bam (in (header num-refs ref-meta) input-file)
    (let* ((hd (subst-sort-order (make-sam-header header) :unsorted))
           (pg (pg-record "readmill-filter"
                          :program-name *readmill-name*
                          :program-version *readmill-version*
                          :previous-program (previous-program hd)
                          :command-line (format nil "~{~a~^ ~}" cmd))))
      (with-bam (out ((with-output-to-string (s)
                        (write-sam-header
                         (add-pg-record hd pg) s)) num-refs ref-meta)
                     output-file :direction :output
                     :if-does-not-exist :create :if-exists :supersede)
        (let* ((counters (mapcar #'make-counting-predicate filters))
               (in (discarding-if (apply #'any counters) in)))
          (loop
             while (has-more-p in)
             do (consume out (next in)))
          (let ((counts (mapcar (lambda (fn)
                                  (rest (multiple-value-list (funcall fn))))
                                counters)))
            (when json-file
              (write-json-file json-file
                               (mapcar (lambda (fn x)
                                         (apply fn x)) descriptors counts)))
            counts))))))

;; Experimental multi-threaded version
(defun pfilter-bam (cmd input-file output-file filters descriptors
                    &optional json-file)
  (with-bam (in (header num-refs ref-meta) input-file)
    (let* ((hd (make-sam-header header))
           (pg (pg-record "readmill-filter"
                          :program-name *readmill-name*
                          :program-version *readmill-version*
                          :previous-program (previous-program hd)
                          :command-line (format nil "~{~a~^ ~}" cmd))))
      (with-bam (out ((with-output-to-string (s)
                        (write-sam-header
                         (add-pg-record hd pg) s)) num-refs ref-meta)
                     output-file :direction :output
                     :if-does-not-exist :create :if-exists :overwrite)
        (let ((counts (batch-filter in out filters
                                    :threads 2 :batch-size 10000)))
          (when json-file
            (write-json-file json-file
                             (mapcar (lambda (fn x)
                                       (apply fn x)) descriptors counts)))
          counts)))))

;; Note that the filter predicates cause items to be removed when they
;; return T
(defun describe-filter-result (num-calls num-filtered name
                               &optional description)
  "Returns an alist mapping the keys NAME DESCRIPTION NUM-IN
NUM-PASSED and NUM-FAILED to the appropriate argument value."
  (pairlis '(name description num-in num-passed num-failed)
           (list name description num-calls (- num-calls num-filtered)
                 num-filtered)))
