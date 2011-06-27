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

(in-package :cl-user)

;; Assumes a recent ASDF with configuration API
(require 'asdf)

;; ;; Ignore any installed libraries; this build should be self-contained
(asdf:clear-configuration)

;; Set the source registry to look in the lib directory only
(asdf:initialize-source-registry
 `(:source-registry
   ,@(mapcar (lambda (path)
               (list :directory (merge-pathnames path)))
             '("lib/deoxybyte-bundle/alexandria/"
               "lib/deoxybyte-bundle/babel/"
               "lib/deoxybyte-bundle/cffi/"
               "lib/deoxybyte-bundle/cl-fad-0.6.4/"
               "lib/deoxybyte-bundle/deoxybyte-gzip/"
               "lib/deoxybyte-bundle/deoxybyte-io/"
               "lib/deoxybyte-bundle/deoxybyte-run/"
               "lib/deoxybyte-bundle/deoxybyte-systems/"
               "lib/deoxybyte-bundle/deoxybyte-unix/"
               "lib/deoxybyte-bundle/deoxybyte-utilities/"
               "lib/deoxybyte-bundle/getopt/"
               "lib/deoxybyte-bundle/trivial-features/"
               "lib/bordeaux-threads/"
               "lib/cl-json/"
               "lib/cl-sam/"
               "lib/eager-future/"
               "."))
   :ignore-inherited-configuration))

;; Set the output translation to put fasl files in the build directory
(asdf:initialize-output-translations
 (list :output-translations (list t (merge-pathnames "build/**/*.*"))
       :ignore-inherited-configuration))

(asdf:load-system :deoxybyte-systems)
(asdf:load-system :readmill)
(asdf:clear-configuration)

(sb-ext:save-lisp-and-die "readmill"
                          :executable t
                          :save-runtime-options t
                          :toplevel #'readmill:readmill-cli)
