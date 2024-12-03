pure function convert_cstr(cstr) result(fstr)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cstr(*)
   character(len=256) :: fstr
   integer :: i

   fstr = ""

   do i = 1, 256
      if (cstr(i) == c_null_char) then
         exit
      end if
      fstr(i:i) = cstr(i)
   end do

   return
end function convert_cstr

subroutine convert_fstr(fstr, cstr)
   use iso_c_binding
   implicit none

   character(*) :: fstr
   character(len=1, kind=c_char), intent(inout) :: cstr(*)

   integer :: i, n

   n = len_trim(fstr)
   do i = 1, n
      cstr(i) = fstr(i:i)
   end do
   cstr(n + 1) = c_null_char

end subroutine convert_fstr

function fsync(fd) bind(c, name="fsync_")
   use iso_c_binding
   integer(c_int), value :: fd
   integer(c_int) :: fsync
end function fsync

subroutine fc_initlagrit(mode, log_file, batch_file) bind(c)
   use iso_c_binding
   implicit none

   integer(c_int32_t), intent(in) :: mode
   character(len=1, kind=c_char), intent(in) :: log_file(*), batch_file(*)

   character(len=8) :: mode_str

   ! general function
   character(len=256) :: convert_cstr

   character(len=256) :: f_log_file, f_batch_file

   if (mode == 0) then
      mode_str = 'silent'
   else if (mode == 1) then
      mode_str = 'noisy'
   end if

   f_log_file = convert_cstr(log_file)
   f_batch_file = convert_cstr(batch_file)
   call initlagrit(mode_str, f_log_file, f_batch_file)

end subroutine fc_initlagrit

subroutine fc_fflush_and_sync(file) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: file(*)

   ! general function
   character(len=256) :: convert_cstr
   integer(c_int) :: fsync

   integer :: fid, fsync_ret

   inquire (file=trim(convert_cstr(file)), number=fid)
   flush (fid)
   fsync_ret = fsync(fnum(fid))

end subroutine fc_fflush_and_sync

subroutine fc_fclose(file) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: file(*)

   ! general function
   character(len=256) :: convert_cstr

   integer :: fid

   inquire (file=trim(convert_cstr(file)), number=fid)
   close (fid)

end subroutine fc_fclose

subroutine fc_dotask(cmd, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cmd
   integer(c_int32_t), intent(inout) :: status

   ! general function
   character(len=256) :: convert_cstr

   call dotask(trim(convert_cstr(cmd))//'; finish', status)

end subroutine fc_dotask

subroutine fc_attr_len(aname, pname, arr_len, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: aname(*), pname(*)
   integer(c_int32_t), intent(inout) :: arr_len, status

   pointer(ptr, buffer)
   integer(c_int8_t) :: buffer(*)

   ! general function
   character(len=256) :: convert_cstr

   call mmfindbk(convert_cstr(aname), convert_cstr(pname), ptr, arr_len, status)

end subroutine fc_attr_len

subroutine fc_mmrelprt(pname, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: pname(*)
   integer(c_int32_t), intent(inout) :: status

   ! general function
   character(len=256) :: convert_cstr

   call mmrelprt(convert_cstr(pname), status)

end subroutine fc_mmrelprt

subroutine fc_mmfindbk(byte_len, cell_len, aname, pname, aptr, arr_len, status) bind(c)
   use iso_c_binding
   implicit none

   integer(c_size_t), intent(in) :: byte_len, cell_len
   character(len=1, kind=c_char), intent(in) :: aname(*), pname(*)
   integer(c_int8_t), intent(inout) :: aptr(*)
   integer(c_int32_t), intent(inout) :: arr_len, status

   pointer(ptr, buffer)
   integer(c_int8_t) :: buffer(*)
   integer :: i

   ! general function
   character(len=256) :: convert_cstr

   call mmfindbk(convert_cstr(aname), convert_cstr(pname), ptr, arr_len, status)

   do i = 1, int(byte_len)*int(cell_len)*int(arr_len)
      aptr(i) = buffer(i)
   end do

end subroutine fc_mmfindbk

subroutine fc_cmo_get_index(cmo_name, idx, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cmo_name(*)
   integer(c_int32_t), intent(inout) :: idx, status

   ! general function
   character(len=256) :: convert_cstr

   call cmo_get_index(convert_cstr(cmo_name), idx, status)

end subroutine fc_cmo_get_index

subroutine fc_cmo_get_name(cmo_name, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(inout) :: cmo_name(*)
   integer(c_int32_t), intent(inout) :: status

   character(len=32) :: tmp

   call cmo_get_name(tmp, status)
   call convert_fstr(tmp, cmo_name)

end subroutine fc_cmo_get_name

subroutine fc_cmo_get_mesh_type(cmo_name, mesh_type, imesh_type, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cmo_name(*)
   character(len=1, kind=c_char), intent(inout) :: mesh_type(*)
   integer(c_int32_t), intent(inout) :: imesh_type, status

   ! general function
   character(len=256) :: convert_cstr

   character(len=32) :: tmp

   call cmo_get_mesh_type(convert_cstr(cmo_name), tmp, imesh_type, status)
   call convert_fstr(tmp, mesh_type)

end subroutine fc_cmo_get_mesh_type
