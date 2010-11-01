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

(defun quality-filter-bam (cmd input-file output-file quality-threshold
                           read-start read-end)
  (with-bam (in (header num-refs ref-meta) input-file)
    (with-bam (out (header num-refs ref-meta) output-file :direction :output
                   :if-does-not-exist :create :if-exists :overwrite)
      (let* ((qfilter (make-counting-predicate
                       (make-quality-p quality-threshold
                                       :read-start read-start
                                       :read-end read-end)))
             (in (discarding-if qfilter in)))
        (loop
           while (has-more-p in)
           do (consume out (next in)))
        (funcall qfilter nil)))))
