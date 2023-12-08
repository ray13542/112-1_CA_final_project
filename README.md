# CA_final_project
Spec: 
    ◆All inputs are synchronized with the negative edge clock.

    ◆All outputs should be synchronized at clock rising edge.

    ◆You should reset all your outputs when i_rst_n is low.

Todo:
    *Submodule
    1. ALU
        (i: 2 data o: 1 data)
        
    2. Control unit
        (i: opcode o: control signals)
        input opcode to determine which control signal need to set
    3. Adder for branch

    4. PC

    5. Muldiv

    6. Memory
        stall = 1 such that next_PC = PC
        rdata => Reg_W_data 
        o_DMEM_cen = MemRead | MemWrite
        o_DMEM_wen = MemWrite & (~MemRead)
        o_DMEM_addr = ALU_result
        o_DMEM_wdata = reg_rdata_2



    7. Write Back

    8. Imediate Generation

    9.ALU control unit
        input: ALUop(2bit) func3, func7
        output: ALUctrl
    10. jal : next_pc = pc + imm
        jalr : next_pc = pc + result of alu
    
    *always block
    1. Intruction fetch
        Set high o_IMEM_cen to load instruction
        pass PC address
        return i_IMEM_data(instruction)

    2. PC-FSM
        PC+4 or PC(stall for MUL) or Branch
    
    

Q:
    load store有沒有要signed extention? A:有 把Instruction某部分擷取成32bit

    Mux怎麼放? A: 在main module裡，assign wire XXX = (control signal) ? : ;

    MUL還需不需要傳ALUctrl? 還是直接輸入i_valid代表要乘法就好?

    BEQ、BNE、BLT、BGE在ALU要做甚麼?
    2 signal (zero lessthan) 
    BEQ (when zero = 1)、BNE(when zero = 0)、BLT(when lessthan = 1)、BGE(when lessthan = 0)
