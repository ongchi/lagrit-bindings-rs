! Convert C strings to Fortran strings
pure function fstr(cstr)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cstr(*)

   ! Most of strings in LaGriT are 32 characters long
   character(len=32) :: fstr
   integer :: i

   fstr = ""

   do i = 1, 32
      if (cstr(i) == c_null_char) then
         exit
      end if
      fstr(i:i) = cstr(i)
   end do

   return
end function fstr

! Get the length of a C string
subroutine c_str_len(cstr, len)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cstr(*)
   integer, intent(out) :: len

   len = 0
   do
      if (cstr(len + 1) == c_null_char) then
         exit
      end if
      len = len + 1
   end do

end subroutine c_str_len

! Convert C strings to Fortran strings
subroutine convert_cstr(cstr, fstr)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cstr(*)
   character(*), intent(out) :: fstr
   integer :: i

   do i = 1, len(fstr)
      fstr(i:i) = cstr(i)
   end do

end subroutine convert_cstr

! Convert Fortran strings to C strings
subroutine convert_fstr(fstr, cstr)
   use iso_c_binding
   implicit none

   character(*), intent(in) :: fstr
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

   ! general function
   character(len=32) :: fstr

   character(len=8) :: mode_str

   if (mode == 0) then
      mode_str = 'silent'
   else if (mode == 1) then
      mode_str = 'noisy'
   end if

   call initlagrit(mode_str, trim(fstr(log_file)), trim(fstr(batch_file)))

end subroutine fc_initlagrit

subroutine fc_fflush_and_sync(file) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: file(*)

   ! general function
   character(len=32) :: fstr
   integer(c_int) :: fsync

   integer :: fid, fsync_ret

   inquire (file=trim(fstr(file)), number=fid)
   flush (fid)
   fsync_ret = fsync(fnum(fid))

end subroutine fc_fflush_and_sync

subroutine fc_fclose(file) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: file(*)

   ! general function
   character(len=32) :: fstr

   integer :: fid

   inquire (file=trim(fstr(file)), number=fid)
   close (fid)

end subroutine fc_fclose

subroutine fc_dotask(cmd, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: cmd
   integer(c_int32_t), intent(inout) :: status

   integer :: str_len
   character(len=:), allocatable :: f_cmd

   ! Because of the length of the cmd string varies.
   ! We do not call fstr(cmd) directly.
   call c_str_len(cmd, str_len)
   allocate (character(str_len) :: f_cmd)
   call convert_cstr(cmd, f_cmd)

   call dotask(f_cmd//";finish", status)

end subroutine fc_dotask

subroutine fc_attr_len(aname, pname, arr_len, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: aname(*), pname(*)
   integer(c_int32_t), intent(inout) :: arr_len, status

   ! general function
   character(len=32) :: fstr

   pointer(ptr, buffer)
   integer(c_int8_t) :: buffer(*)

   call mmfindbk(fstr(aname), fstr(pname), ptr, arr_len, status)

end subroutine fc_attr_len

subroutine fc_mmrelprt(pname, status) bind(c)
   use iso_c_binding
   implicit none

   character(len=1, kind=c_char), intent(in) :: pname(*)
   integer(c_int32_t), intent(inout) :: status

   ! general function
   character(len=32) :: fstr

   call mmrelprt(fstr(pname), status)

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
   character(len=32) :: fstr

   call mmfindbk(fstr(aname), fstr(pname), ptr, arr_len, status)

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
   character(len=32) :: fstr

   call cmo_get_index(fstr(cmo_name), idx, status)

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
   character(len=32) :: fstr

   character(len=32) :: tmp

   call cmo_get_mesh_type(fstr(cmo_name), tmp, imesh_type, status)
   call convert_fstr(tmp, mesh_type)

end subroutine fc_cmo_get_mesh_type

subroutine fc_control_command_lg(status) bind(c)
   use iso_c_binding
   implicit none

   integer(c_int32_t), intent(inout) :: status

   call control_command_lg(status)

end subroutine fc_control_command_lg

subroutine fc_set_iattr(attr, cmo_name, data, data_len, status) bind(c)
   ! This subroutine does not check the type of the attribute when writing to buffer.
   use iso_c_binding
   use iso_fortran_env
   implicit none

   character(len=1, kind=c_char), intent(in) :: attr(*), cmo_name(*)
   integer(c_int64_t), intent(in) :: data(*), data_len
   integer(c_int32_t), intent(inout) :: status

   ! general function
   character(len=32) :: fstr

   integer(int64) :: i, attr_len, itype, buffer(*)
   pointer(ipattr, buffer)

   call cmo_get_info(fstr(attr), fstr(cmo_name), ipattr, attr_len, itype, status)

   do i = 1, data_len
      buffer(i) = data(i)
   end do

end subroutine fc_set_iattr

subroutine fc_set_fattr(attr, cmo_name, data, data_len, status) bind(c)
   ! This subroutine does not check the type of the attribute when writing to buffer.
   use iso_c_binding
   use iso_fortran_env
   implicit none

   character(len=1, kind=c_char), intent(in) :: attr(*), cmo_name(*)
   real(c_double), intent(in) :: data(*)
   integer(c_int64_t), intent(in) :: data_len
   integer(c_int32_t), intent(inout) :: status

   ! general function
   character(len=32) :: fstr

   integer(int64) :: i, attr_len, itype
   real(real64) :: buffer(*)
   pointer(ipattr, buffer)

   call cmo_get_info(fstr(attr), fstr(cmo_name), ipattr, attr_len, itype, status)

   do i = 1, data_len
      buffer(i) = data(i)
   end do

end subroutine fc_set_fattr
