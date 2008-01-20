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

    dword() {
        opcodes=( ${opcodes[@]} "oc_$1" )
    }

    qword() {
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

    void() {
        if [[ -n ${2} ]] ; then
            echo "void ${2}(const Instruction &);"
        fi
    }

    dword() {
        if [[ -n ${2} ]] ; then
            echo "unsigned ${2}(const Instruction &);"
        fi
    }

    qword() {
        if [[ -n ${2} ]] ; then
            echo "unsigned long long ${2}(const Instruction &);"
        fi
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
        if [[ -z ${2} ]] ; then
            echo "                    operation(operations::${1}, instructions[instruction_index]);"
        else
            echo "                    ${2}(instructions[instruction_index]);"
        fi
        echo "                    debug_leave();"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    dword() {
        echo "                case oc_${1}:"
        echo "                    debug_enter();"
        if [[ -z ${2} ]] ; then
            echo "                    retval = operation(operations::${1}, instructions[instruction_index]);"
        else
            echo "                    retval = ${2}(instructions[instruction_index]);"
        fi
        echo "                    debug_leave();"
        echo "                    spu_write_out_intr_mbox(km_result_dword);"
        echo "                    spu_write_out_mbox(retval & 0xFFFFFFFF);"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    return_qword() {
        echo "                case oc_${1}:"
        echo "                    debug_enter();"
        if [[ -z ${2} ]] ; then
            echo "                    retval = operation(operations::${1}, instructions[instruction_index]);"
        else
            echo "                    retval = ${2}(instructions[instruction_index]);"
        fi
        echo "                    debug_leave();"
        echo "                    spu_write_out_intr_mbox(km_result_qword);"
        echo "                    spu_write_out_mbox(retval & 0xFFFFFFFF);"
        echo "                    spu_write_out_mbox(retval >> 32);"
        echo "                    spu_write_out_intr_mbox(km_instruction_finished);"
        echo "                    break;"
        echo
    }

    source ${sksource} > ${skbody}
)
