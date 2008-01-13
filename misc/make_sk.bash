#!/bin/bash
# vim: set sw=4 sts=4 et tw=0 :

sksource=${1}
skbody=${1%.sk}.body
skfunctions=${1%.sk}.functions
skopcodes=${1%.sk}.opcodes

(
    opcodes=( )
    inline_function() {
        :
    }

    halt() {
        opcodes=( ${opcodes[@]} "oc_halt" )
    }

    void() {
        opcodes=( ${opcodes[@]} "oc_$1" )
    }

    return_void() {
        opcodes=( ${opcodes[@]} "oc_$1" )
    }

    dword() {
        opcodes=( ${opcodes[@]} "oc_$1" )
    }

    return_dword() {
        opcodes=( ${opcodes[@]} "oc_$1" )
    }

    return_qword() {
        opcodes=( ${opcodes[@]} "oc_$1" )
    }

    source ${sksource}

    echo -e "${opcodes[@]}" > ${skopcodes}

    sed -i \
        -e 's/ /,\n    /g' \
        -e 's/^/    /' \
        ${skopcodes}
)

(
    inline_function() {
        echo "inline ${1} ${2}(Instruction &)"
        echo "{"
        cat
        echo "}"
        echo
    }

    halt() {
        :
    }

    return_void() {
        echo "void ${2}(const Instruction &);"
    }

    void() {
        :
    }

    return_dword() {
        echo "unsigned ${2}(const Instruction &);"
    }

    dword() {
        :
    }

    return_qword() {
        echo "unsigned long long ${2}(const Instruction &);"
    }

    source ${sksource} > ${skfunctions}
)

(
    inline_function() {
        :
    }

    halt() {
        echo "                case oc_halt:"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    spu_stop(0);"
        echo "                    break;"
        echo
    }

    void() {
        echo "                case oc_${1}:"
        echo "                    debug_enter();"
        echo "                    operation(operations::${1}, instructions[instruction_index]);"
        echo "                    debug_leave();"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    return_void() {
        echo ">>> return_void is deprecated syntax! <<<" >&2
        echo ">>>     Please use 'void' instead!    <<<" >&2
        echo "                case oc_${1}:"
        echo "                    #ifdef DEBUG"
        echo "                    spu_write_out_intr_mbox(km_debug_enter);"
        echo "                    #endif"
        echo "                    ${2}(instructions[instruction_index]);"
        echo "                    #ifdef DEBUG"
        echo "                    spu_write_out_intr_mbox(km_debug_leave);"
        echo "                    #endif"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    dword() {
        echo "                case oc_${1}:"
        echo "                    debug_enter();"
        echo "                    retval = operation(operations::${1}, instructions[instruction_index]);"
        echo "                    debug_leave();"
        echo "                    spu_write_out_intr_mbox(km_result_dword);"
        echo "                    spu_write_out_mbox(retval & 0xFFFFFFFF);"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    return_dword() {
        echo ">>> return_dword is deprecated syntax! <<<" >&2
        echo ">>>     Please use 'dword' instead!    <<<" >&2
        echo "                case oc_${1}:"
        echo "                    #ifdef DEBUG"
        echo "                    spu_write_out_intr_mbox(km_debug_enter);"
        echo "                    #endif"
        echo "                    retval = ${2}(instructions[instruction_index]);"
        echo "                    #ifdef DEBUG"
        echo "                    spu_write_out_intr_mbox(km_debug_leave);"
        echo "                    #endif"
        echo "                    spu_write_out_intr_mbox(km_result_dword);"
        echo "                    spu_write_out_mbox(retval & 0xFFFFFFFF);"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    return_qword() {
        echo "                case oc_${1}:"
        echo "                    #ifdef DEBUG"
        echo "                    spu_write_out_intr_mbox(km_debug_enter);"
        echo "                    #endif"
        echo "                    retval = ${2}(instructions[instruction_index]);"
        echo "                    #ifdef DEBUG"
        echo "                    spu_write_out_intr_mbox(km_debug_leave);"
        echo "                    #endif"
        echo "                    spu_write_out_intr_mbox(km_result_dword);"
        echo "                    spu_write_out_mbox(retval & 0xFFFFFFFF);"
        echo "                    spu_write_out_mbox(retval >> 32);"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    source ${sksource} > ${skbody}
)
